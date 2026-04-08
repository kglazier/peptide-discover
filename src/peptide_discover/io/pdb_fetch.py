"""RCSB PDB file fetching."""

from pathlib import Path

import requests

RCSB_API = "https://files.rcsb.org/download"


def fetch_pdb(pdb_id: str, output_dir: Path) -> Path:
    """Download a PDB file from RCSB."""
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{pdb_id.lower()}.pdb"

    if output_path.exists():
        return output_path

    url = f"{RCSB_API}/{pdb_id.upper()}.pdb"
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()

    output_path.write_text(resp.text)
    return output_path
