"""Stage 3: Binding affinity prediction via AutoDock Vina."""

import json
import logging
import subprocess
import tempfile
from pathlib import Path

from peptide_discover.models.peptide import PeptideCandidate
from peptide_discover.models.scores import BindingResult
from peptide_discover.models.target import TargetProtein

logger = logging.getLogger(__name__)

VINA_BINARY = Path(__file__).parents[3] / "bin" / "vina.exe"


def predict_binding(
    target: TargetProtein,
    candidates: list[PeptideCandidate],
    top_k: int = 50,
    exhaustiveness: int = 8,
) -> list[BindingResult]:
    """Predict binding affinity using AutoDock Vina.

    Docks each peptide candidate against the target protein structure
    and returns top_k candidates sorted by binding energy (kcal/mol).

    Args:
        target: Target protein with structure file.
        candidates: Peptide candidates to dock.
        top_k: Keep top K candidates by binding score.
        exhaustiveness: Vina exhaustiveness (higher = more thorough, slower).
    """
    if not target.structure_path or not target.structure_path.exists():
        raise ValueError(
            f"Target '{target.identifier}' has no structure file. "
            "Binding prediction requires a 3D structure."
        )

    vina_path = _find_vina()

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Prepare receptor PDBQT
        logger.info("Preparing receptor...")
        receptor_pdbqt = _prepare_receptor(target.structure_path, tmpdir)

        # Get receptor center and box size
        center, box_size = _compute_binding_box(target.structure_path)
        logger.info(
            "Search box: center=(%.1f, %.1f, %.1f), size=(%.1f, %.1f, %.1f)",
            *center, *box_size
        )

        # Dock each candidate
        results = []
        for i, candidate in enumerate(candidates):
            logger.info(
                "Docking candidate %d/%d: %s",
                i + 1, len(candidates), candidate.sequence,
            )
            try:
                result = _dock_peptide(
                    candidate=candidate,
                    receptor_pdbqt=receptor_pdbqt,
                    target_id=target.identifier,
                    center=center,
                    box_size=box_size,
                    vina_path=vina_path,
                    exhaustiveness=exhaustiveness,
                    work_dir=tmpdir,
                )
                results.append(result)
            except Exception as e:
                logger.warning("Failed to dock %s: %s", candidate.sequence, e)
                results.append(
                    BindingResult(
                        peptide_sequence=candidate.sequence,
                        target_id=target.identifier,
                        affinity_score=0.0,
                        confidence=0.0,
                    )
                )

    # Sort by affinity (more negative = better binding)
    results.sort(key=lambda r: r.affinity_score)

    # Keep top_k
    results = results[:top_k]

    logger.info(
        "Binding prediction complete. Top score: %.2f kcal/mol",
        results[0].affinity_score if results else 0,
    )
    return results


def _find_vina() -> Path:
    """Locate the Vina binary."""
    if VINA_BINARY.exists():
        return VINA_BINARY

    # Check PATH
    import shutil
    vina = shutil.which("vina") or shutil.which("vina.exe")
    if vina:
        return Path(vina)

    raise FileNotFoundError(
        f"AutoDock Vina not found at {VINA_BINARY} or on PATH. "
        "Download from: https://github.com/ccsb-scripps/AutoDock-Vina/releases"
    )


def _prepare_receptor(pdb_path: Path, work_dir: Path) -> Path:
    """Convert receptor PDB to PDBQT format using OpenBabel.

    OpenBabel handles hydrogen addition, charge assignment, and
    proper AutoDock atom typing — the standard approach for Vina.
    """
    from openbabel import openbabel

    output_path = work_dir / "receptor.pdbqt"

    conv = openbabel.OBConversion()
    conv.SetInFormat("pdb")
    conv.SetOutFormat("pdbqt")
    # Add hydrogens, remove waters
    conv.AddOption("r", conv.OUTOPTIONS)  # rigid (no rotatable bonds for receptor)
    conv.AddOption("x", conv.OUTOPTIONS)  # remove non-polar hydrogens from output
    conv.AddOption("p", conv.OUTOPTIONS)  # add hydrogens at pH 7.4

    mol = openbabel.OBMol()
    conv.ReadFile(mol, str(pdb_path))

    # Remove water molecules
    waters = []
    for res in openbabel.OBResidueIter(mol):
        if res.GetName().strip() in ("HOH", "WAT"):
            waters.append(res)
    for w in waters:
        mol.DeleteResidue(w)

    mol.AddHydrogens(False, True, 7.4)  # polar only, pH 7.4

    conv.WriteFile(mol, str(output_path))

    # Count atoms in output
    atom_count = sum(1 for line in open(output_path) if line.startswith("ATOM"))
    logger.info("Receptor PDBQT: %d atoms", atom_count)
    return output_path


def _compute_binding_box(pdb_path: Path) -> tuple[tuple, tuple]:
    """Compute a search box that covers the entire protein surface.

    Returns (center_xyz, size_xyz) in Angstroms.
    """
    import numpy as np

    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])

    coords = np.array(coords)
    center = coords.mean(axis=0)
    extent = coords.max(axis=0) - coords.min(axis=0)

    # Add padding for peptide binding at surface
    box_size = extent + 10.0  # 5A padding on each side

    # Large boxes massively slow docking. Cap at 50A — covers most
    # binding sites while keeping search feasible on consumer hardware.
    # For full-protein screening, multiple smaller boxes would be better,
    # but for MVP this is a reasonable tradeoff.
    box_size = np.clip(box_size, 20.0, 50.0)

    return tuple(center), tuple(box_size)


def _sequence_to_pdbqt(sequence: str, work_dir: Path) -> Path:
    """Generate 3D peptide structure from sequence and convert to PDBQT."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Build peptide SMILES from amino acid sequence
    smiles = _peptide_to_smiles(sequence)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Failed to parse peptide SMILES for sequence: {sequence}")

    # Add hydrogens and generate 3D conformation
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        # Fallback: try with random coordinates
        params.useRandomCoords = True
        result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            raise ValueError(f"Failed to generate 3D conformation for: {sequence}")

    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)

    # Use meeko to convert to PDBQT
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    mol_setup = mol_setups[0]

    pdbqt_string, is_ok, err_msg = PDBQTWriterLegacy.write_string(mol_setup)
    if not is_ok:
        raise ValueError(f"PDBQT conversion failed for {sequence}: {err_msg}")

    output_path = work_dir / f"peptide_{sequence[:8]}.pdbqt"
    with open(output_path, "w") as f:
        f.write(pdbqt_string)

    return output_path


def _peptide_to_smiles(sequence: str) -> str:
    """Convert amino acid sequence to SMILES string."""
    from rdkit import Chem

    # Use RDKit's built-in peptide builder
    mol = Chem.MolFromSequence(sequence)
    if mol is None:
        raise ValueError(f"Cannot convert sequence to molecule: {sequence}")
    return Chem.MolToSmiles(mol)


def _dock_peptide(
    candidate: PeptideCandidate,
    receptor_pdbqt: Path,
    target_id: str,
    center: tuple,
    box_size: tuple,
    vina_path: Path,
    exhaustiveness: int,
    work_dir: Path,
) -> BindingResult:
    """Dock a single peptide against the receptor."""
    # Generate peptide PDBQT
    ligand_pdbqt = _sequence_to_pdbqt(candidate.sequence, work_dir)
    output_pdbqt = work_dir / f"out_{candidate.sequence[:8]}.pdbqt"

    # Run Vina
    cmd = [
        str(vina_path),
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--out", str(output_pdbqt),
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{box_size[0]:.3f}",
        "--size_y", f"{box_size[1]:.3f}",
        "--size_z", f"{box_size[2]:.3f}",
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", "1",
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=600,
    )

    if result.returncode != 0:
        raise RuntimeError(f"Vina failed: {result.stderr}")

    # Parse best binding energy from output
    affinity = _parse_vina_output(result.stdout)

    return BindingResult(
        peptide_sequence=candidate.sequence,
        target_id=target_id,
        affinity_score=affinity,
        confidence=min(1.0, abs(affinity) / 15.0),  # rough confidence normalization
        complex_pdb_path=str(output_pdbqt) if output_pdbqt.exists() else None,
    )


def _parse_vina_output(stdout: str) -> float:
    """Parse the best binding energy from Vina stdout."""
    # Look for the results table after the header line
    in_results = False
    for line in stdout.split("\n"):
        stripped = line.strip()
        if "-----+----" in stripped:
            in_results = True
            continue
        if in_results and stripped:
            # Format: "   1       -8.093          0          0"
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    mode = int(parts[0])
                    affinity = float(parts[1])
                    return affinity
                except ValueError:
                    continue

    raise RuntimeError("Could not parse binding energy from Vina output")
