"""Stage 1: Target protein preparation."""

from pathlib import Path

from peptide_discover.models.target import TargetProtein


def resolve_target(identifier: str, output_dir: Path = Path("data/targets")) -> TargetProtein:
    """Resolve a target protein from UniProt ID, PDB ID, or file path.

    Fetches sequence and structure, stores locally.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Detect identifier type
    path = Path(identifier)
    if path.exists():
        return _from_file(path)

    # UniProt IDs are typically 6-10 alphanumeric chars
    if len(identifier) >= 6 and identifier[0] in "OPQA":
        return _from_uniprot(identifier, output_dir)

    # PDB IDs are 4 chars
    if len(identifier) == 4:
        return _from_pdb(identifier, output_dir)

    raise ValueError(
        f"Cannot resolve '{identifier}'. Provide a UniProt ID, PDB ID, or file path."
    )


def _from_file(path: Path) -> TargetProtein:
    """Load target from a local PDB or FASTA file."""
    from peptide_discover.io.formats import read_fasta, read_pdb_sequence

    if path.suffix in (".fasta", ".fa", ".faa"):
        sequence = read_fasta(path)
        return TargetProtein(identifier=path.name, sequence=sequence, structure_path=None)
    elif path.suffix in (".pdb", ".cif"):
        sequence = read_pdb_sequence(path)
        return TargetProtein(
            identifier=path.name, sequence=sequence, structure_path=path
        )
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")


def _from_uniprot(uniprot_id: str, output_dir: Path) -> TargetProtein:
    """Fetch target from UniProt + AlphaFold DB."""
    from peptide_discover.io.uniprot import fetch_sequence, fetch_alphafold_structure

    sequence = fetch_sequence(uniprot_id)
    structure_path = fetch_alphafold_structure(uniprot_id, output_dir)

    return TargetProtein(
        identifier=uniprot_id,
        uniprot_id=uniprot_id,
        sequence=sequence,
        structure_path=structure_path,
    )


def _from_pdb(pdb_id: str, output_dir: Path) -> TargetProtein:
    """Fetch target from RCSB PDB."""
    from peptide_discover.io.pdb_fetch import fetch_pdb

    structure_path = fetch_pdb(pdb_id, output_dir)

    return TargetProtein(
        identifier=pdb_id,
        pdb_id=pdb_id,
        structure_path=structure_path,
    )
