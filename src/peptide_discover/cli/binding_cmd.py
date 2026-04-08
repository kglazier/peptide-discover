"""CLI commands for binding prediction."""

from pathlib import Path

import typer

app = typer.Typer()


@app.command("predict")
def predict(
    candidates_file: Path = typer.Argument(help="Path to generated candidates (CSV/JSON)."),
    target: str = typer.Option(..., "--target", "-t", help="Target protein identifier."),
    top_k: int = typer.Option(50, "--top-k", "-k"),
) -> None:
    """Predict binding affinity for peptide candidates using Boltz-2."""
    from peptide_discover.stages.binding import predict_binding

    typer.echo(f"Binding prediction not yet implemented. Target: {target}, file: {candidates_file}")
