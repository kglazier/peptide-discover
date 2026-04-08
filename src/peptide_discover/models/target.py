"""Target protein data models."""

from pathlib import Path

from pydantic import BaseModel


class BindingSite(BaseModel):
    """A binding site on the target protein."""

    chain: str = "A"
    residue_indices: list[int] = []
    description: str = ""


class TargetProtein(BaseModel):
    """A target protein for peptide design."""

    identifier: str  # UniProt ID, PDB ID, or file path
    name: str = ""
    sequence: str = ""
    structure_path: Path | None = None
    uniprot_id: str | None = None
    pdb_id: str | None = None
    binding_sites: list[BindingSite] = []
