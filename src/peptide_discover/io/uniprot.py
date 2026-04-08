"""UniProt REST API client."""

from pathlib import Path

import requests

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api"


def fetch_sequence(uniprot_id: str) -> str:
    """Fetch protein sequence from UniProt."""
    url = f"{UNIPROT_API}/{uniprot_id}.fasta"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()

    lines = resp.text.strip().split("\n")
    return "".join(line for line in lines if not line.startswith(">"))


def fetch_alphafold_structure(uniprot_id: str, output_dir: Path) -> Path | None:
    """Download predicted structure from AlphaFold DB.

    Uses the AlphaFold API to resolve the current model version.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check for any existing cached version
    existing = list(output_dir.glob(f"AF-{uniprot_id}-F1-model_v*.pdb"))
    if existing:
        return existing[0]

    # Query API for the correct PDB URL
    api_url = f"{ALPHAFOLD_API}/prediction/{uniprot_id}"
    api_resp = requests.get(api_url, timeout=30)

    if api_resp.status_code == 404:
        return None

    api_resp.raise_for_status()
    entries = api_resp.json()
    if not entries:
        return None

    pdb_url = entries[0].get("pdbUrl")
    if not pdb_url:
        return None

    # Extract version from URL for filename
    version = entries[0].get("latestVersion", "unknown")
    output_path = output_dir / f"AF-{uniprot_id}-F1-model_v{version}.pdb"

    resp = requests.get(pdb_url, timeout=60)
    resp.raise_for_status()
    output_path.write_text(resp.text)
    return output_path
