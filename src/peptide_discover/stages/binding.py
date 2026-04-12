"""Stage 3: Binding affinity prediction.

Supports two docking engines:
- vina: AutoDock Vina (fast, best for 8-mers, native Windows)
- adcp: AutoDock CrankPep (handles 15+ mers, folds and docks, runs via WSL)
"""

import logging
import re
import subprocess
import tempfile
from pathlib import Path

from peptide_discover.models.peptide import PeptideCandidate
from peptide_discover.models.scores import BindingResult
from peptide_discover.models.target import TargetProtein

logger = logging.getLogger(__name__)

VINA_BINARY = Path(__file__).parents[3] / "bin" / "vina.exe"
ADCP_CACHE_DIR = Path(__file__).parents[3] / "data" / "adcp_work"

# WSL preamble to activate ADCP environment
WSL_ADCP_PREAMBLE = (
    "export MAMBA_EXE=/home/kglazier/.local/bin/micromamba && "
    "export MAMBA_ROOT_PREFIX=/home/kglazier/micromamba && "
    'eval "$($MAMBA_EXE shell hook --shell bash --root-prefix $MAMBA_ROOT_PREFIX 2>/dev/null)" && '
    "micromamba activate adcpsuite"
)


def predict_binding(
    target: TargetProtein,
    candidates: list[PeptideCandidate],
    top_k: int = 50,
    exhaustiveness: int = 8,
    engine: str = "vina",
    adcp_runs: int = 10,
    adcp_steps: int = 2500000,
) -> list[BindingResult]:
    """Predict binding affinity for peptide candidates.

    Args:
        target: Target protein with structure file.
        candidates: Peptide candidates to dock.
        top_k: Keep top K candidates by binding score.
        exhaustiveness: Vina exhaustiveness (ignored for ADCP).
        engine: "vina" or "adcp".
        adcp_runs: Number of ADCP MC runs per peptide (more = better sampling).
        adcp_steps: MC steps per ADCP run (more = more thorough).
    """
    if not target.structure_path or not target.structure_path.exists():
        raise ValueError(
            f"Target '{target.identifier}' has no structure file. "
            "Binding prediction requires a 3D structure."
        )

    if engine == "adcp":
        results = _predict_binding_adcp(
            target, candidates, adcp_runs=adcp_runs, adcp_steps=adcp_steps,
        )
    elif engine == "vina":
        results = _predict_binding_vina(
            target, candidates, exhaustiveness=exhaustiveness,
        )
    else:
        raise ValueError(f"Unknown docking engine: {engine}. Use 'vina' or 'adcp'.")

    # Sort by affinity (more negative = better binding)
    results.sort(key=lambda r: r.affinity_score)

    # Keep top_k
    results = results[:top_k]

    logger.info(
        "Binding prediction complete. Top score: %.2f kcal/mol",
        results[0].affinity_score if results else 0,
    )
    return results


# =============================================================================
# ADCP (CrankPep) engine — handles long peptides via WSL
# =============================================================================

def _predict_binding_adcp(
    target: TargetProtein,
    candidates: list[PeptideCandidate],
    adcp_runs: int = 10,
    adcp_steps: int = 2500000,
) -> list[BindingResult]:
    """Dock candidates using AutoDock CrankPep via WSL."""
    # Prepare target file (cached per target)
    target_dir = ADCP_CACHE_DIR / target.identifier
    target_dir.mkdir(parents=True, exist_ok=True)
    trg_file = target_dir / f"{target.identifier}_target.trg"

    if not trg_file.exists():
        logger.info("Preparing ADCP target file for %s (this takes ~5-15 min)...", target.identifier)
        _prepare_adcp_target(target.structure_path, target_dir, target.identifier)
    else:
        logger.info("Using cached ADCP target file for %s", target.identifier)

    # Dock each candidate
    results = []
    for i, candidate in enumerate(candidates):
        logger.info(
            "Docking candidate %d/%d: %s (ADCP, %d runs)",
            i + 1, len(candidates), candidate.sequence, adcp_runs,
        )
        try:
            result = _dock_peptide_adcp(
                candidate=candidate,
                trg_file=trg_file,
                target_id=target.identifier,
                target_dir=target_dir,
                adcp_runs=adcp_runs,
                adcp_steps=adcp_steps,
            )
            results.append(result)
            logger.info(
                "  %s: %.2f kcal/mol", candidate.sequence, result.affinity_score,
            )
        except Exception as e:
            logger.warning("Failed to dock %s: %s", candidate.sequence, str(e)[:100])
            results.append(
                BindingResult(
                    peptide_sequence=candidate.sequence,
                    target_id=target.identifier,
                    affinity_score=0.0,
                    confidence=0.0,
                )
            )

    return results


def _win_to_wsl_path(win_path: Path) -> str:
    """Convert a Windows path to WSL path."""
    path_str = str(win_path.resolve()).replace("\\", "/")
    # C:/Users/... -> /mnt/c/Users/...
    if len(path_str) >= 2 and path_str[1] == ":":
        drive = path_str[0].lower()
        return f"/mnt/{drive}{path_str[2:]}"
    return path_str


def _run_wsl_adcp(cmd: str, timeout: int = 1800) -> subprocess.CompletedProcess:
    """Run a command in WSL with ADCP environment activated."""
    full_cmd = f"{WSL_ADCP_PREAMBLE} && {cmd}"
    return subprocess.run(
        ["wsl", "-d", "Ubuntu", "-e", "bash", "-c", full_cmd],
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def _prepare_adcp_target(pdb_path: Path, target_dir: Path, target_name: str) -> None:
    """Prepare ADCP target file (.trg) via WSL.

    Runs: reduce → prepare_receptor → agfr
    """
    wsl_pdb = _win_to_wsl_path(pdb_path)
    wsl_dir = _win_to_wsl_path(target_dir)

    cmd = (
        f"cd {wsl_dir} && "
        f"reduce {wsl_pdb} > receptor_h.pdb 2>/dev/null && "
        f"prepare_receptor -r receptor_h.pdb -o receptor.pdbqt 2>&1 && "
        f"agfr -r receptor.pdbqt -P 10 -o {target_name}_target 2>&1"
    )

    result = _run_wsl_adcp(cmd, timeout=1800)

    trg_file = target_dir / f"{target_name}_target.trg"
    if not trg_file.exists():
        raise RuntimeError(
            f"ADCP target file not created. stdout: {result.stdout[-500:]}\n"
            f"stderr: {result.stderr[-500:]}"
        )

    logger.info("ADCP target file created: %s", trg_file)


def _dock_peptide_adcp(
    candidate: PeptideCandidate,
    trg_file: Path,
    target_id: str,
    target_dir: Path,
    adcp_runs: int,
    adcp_steps: int,
) -> BindingResult:
    """Dock a single peptide using ADCP."""
    wsl_trg = _win_to_wsl_path(trg_file)
    wsl_dir = _win_to_wsl_path(target_dir)

    # ADCP uses lowercase for coil, uppercase for helix
    # Default: all lowercase (coil start) — let ADCP find the structure
    sequence = candidate.sequence.lower()
    job_name = f"dock_{candidate.sequence[:8]}"

    cmd = (
        f"cd {wsl_dir} && "
        f"adcp -T {wsl_trg} -s {sequence} "
        f"-N {adcp_runs} -n {adcp_steps} "
        f"-o {job_name} -w {wsl_dir} -O 2>&1"
    )

    result = _run_wsl_adcp(cmd, timeout=3600)

    if result.returncode != 0:
        raise RuntimeError(f"ADCP failed: {result.stderr[-300:]}")

    # Parse best affinity from output
    affinity = _parse_adcp_output(result.stdout)

    return BindingResult(
        peptide_sequence=candidate.sequence,
        target_id=target_id,
        affinity_score=affinity,
        confidence=min(1.0, abs(affinity) / 15.0),
    )


def _parse_adcp_output(stdout: str) -> float:
    """Parse the best binding affinity from ADCP stdout.

    ADCP output format:
    mode |  affinity  | ref. | clust. | rmsd | energy | best |
         | (kcal/mol) | fnc  |  size  | stdv |  stdv  | run  |
    -----+------------+------+--------+------+--------+------+
       1        -12.4      0.0      11      NA      NA    341
    """
    in_results = False
    for line in stdout.split("\n"):
        stripped = line.strip()
        if "-----+----" in stripped:
            in_results = True
            continue
        if in_results and stripped:
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    int(parts[0])  # mode number
                    return float(parts[1])
                except ValueError:
                    continue

    raise RuntimeError("Could not parse binding energy from ADCP output")


# =============================================================================
# Vina engine — fast, best for short peptides (8-mers)
# =============================================================================

def _predict_binding_vina(
    target: TargetProtein,
    candidates: list[PeptideCandidate],
    exhaustiveness: int = 8,
) -> list[BindingResult]:
    """Dock candidates using AutoDock Vina."""
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
                result = _dock_peptide_vina(
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
                logger.info(
                    "  %s: %.2f kcal/mol", candidate.sequence, result.affinity_score,
                )
            except Exception as e:
                logger.warning("Failed to dock %s: %s", candidate.sequence, str(e)[:100])
                results.append(
                    BindingResult(
                        peptide_sequence=candidate.sequence,
                        target_id=target.identifier,
                        affinity_score=0.0,
                        confidence=0.0,
                    )
                )

    return results


def _find_vina() -> Path:
    """Locate the Vina binary."""
    if VINA_BINARY.exists():
        return VINA_BINARY

    import shutil
    vina = shutil.which("vina") or shutil.which("vina.exe")
    if vina:
        return Path(vina)

    raise FileNotFoundError(
        f"AutoDock Vina not found at {VINA_BINARY} or on PATH. "
        "Download from: https://github.com/ccsb-scripps/AutoDock-Vina/releases"
    )


def _prepare_receptor(pdb_path: Path, work_dir: Path) -> Path:
    """Convert receptor PDB to PDBQT format using OpenBabel."""
    from openbabel import openbabel

    output_path = work_dir / "receptor.pdbqt"

    conv = openbabel.OBConversion()
    conv.SetInFormat("pdb")
    conv.SetOutFormat("pdbqt")
    conv.AddOption("r", conv.OUTOPTIONS)
    conv.AddOption("x", conv.OUTOPTIONS)
    conv.AddOption("p", conv.OUTOPTIONS)

    mol = openbabel.OBMol()
    conv.ReadFile(mol, str(pdb_path))

    waters = []
    for res in openbabel.OBResidueIter(mol):
        if res.GetName().strip() in ("HOH", "WAT"):
            waters.append(res)
    for w in waters:
        mol.DeleteResidue(w)

    mol.AddHydrogens(False, True, 7.4)
    conv.WriteFile(mol, str(output_path))

    atom_count = sum(1 for line in open(output_path) if line.startswith("ATOM"))
    logger.info("Receptor PDBQT: %d atoms", atom_count)
    return output_path


def _compute_binding_box(pdb_path: Path) -> tuple[tuple, tuple]:
    """Compute a search box covering the protein surface."""
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
    box_size = extent + 10.0
    box_size = np.clip(box_size, 20.0, 50.0)

    return tuple(center), tuple(box_size)


def _sequence_to_pdbqt(sequence: str, work_dir: Path) -> Path:
    """Generate 3D peptide structure from sequence and convert to PDBQT."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    smiles = _peptide_to_smiles(sequence)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Failed to parse peptide SMILES for sequence: {sequence}")

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        params.useRandomCoords = True
        result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            raise ValueError(f"Failed to generate 3D conformation for: {sequence}")

    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)

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

    mol = Chem.MolFromSequence(sequence)
    if mol is None:
        raise ValueError(f"Cannot convert sequence to molecule: {sequence}")
    return Chem.MolToSmiles(mol)


def _dock_peptide_vina(
    candidate: PeptideCandidate,
    receptor_pdbqt: Path,
    target_id: str,
    center: tuple,
    box_size: tuple,
    vina_path: Path,
    exhaustiveness: int,
    work_dir: Path,
) -> BindingResult:
    """Dock a single peptide using Vina."""
    ligand_pdbqt = _sequence_to_pdbqt(candidate.sequence, work_dir)
    output_pdbqt = work_dir / f"out_{candidate.sequence[:8]}.pdbqt"

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

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

    if result.returncode != 0:
        raise RuntimeError(f"Vina failed: {result.stderr}")

    affinity = _parse_vina_output(result.stdout)

    return BindingResult(
        peptide_sequence=candidate.sequence,
        target_id=target_id,
        affinity_score=affinity,
        confidence=min(1.0, abs(affinity) / 15.0),
        complex_pdb_path=str(output_pdbqt) if output_pdbqt.exists() else None,
    )


def _parse_vina_output(stdout: str) -> float:
    """Parse the best binding energy from Vina stdout."""
    in_results = False
    for line in stdout.split("\n"):
        stripped = line.strip()
        if "-----+----" in stripped:
            in_results = True
            continue
        if in_results and stripped:
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    int(parts[0])
                    return float(parts[1])
                except ValueError:
                    continue

    raise RuntimeError("Could not parse binding energy from Vina output")
