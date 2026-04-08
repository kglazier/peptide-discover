"""CLI commands for target protein preparation."""

from pathlib import Path

import typer

app = typer.Typer()


@app.command("fetch")
def fetch(
    identifier: str = typer.Argument(help="UniProt ID, PDB ID, or FASTA file path."),
    output_dir: Path = typer.Option(Path("data/targets"), "--output-dir", "-o"),
) -> None:
    """Fetch and prepare a target protein structure."""
    from peptide_discover.stages.target_prep import resolve_target

    target = resolve_target(identifier, output_dir=output_dir)
    typer.echo(f"Target: {target.name or target.identifier}")
    typer.echo(f"Sequence length: {len(target.sequence)}")
    if target.structure_path:
        typer.echo(f"Structure: {target.structure_path}")
