"""CLI commands for ranking and export."""

from pathlib import Path

import typer

app = typer.Typer()


@app.command("export")
def export(
    candidates_file: Path = typer.Argument(help="Path to screened candidates."),
    output: Path = typer.Option(Path("results/ranked.csv"), "--output", "-o"),
    fmt: str = typer.Option("csv", "--format", "-f", help="Output format: csv or json."),
) -> None:
    """Rank candidates and export results."""
    typer.echo(f"Ranking not yet implemented. File: {candidates_file}")
