"""CLI commands for property/safety screening."""

from pathlib import Path

import typer

app = typer.Typer()


@app.command("run")
def run_screening(
    candidates_file: Path = typer.Argument(help="Path to candidates with binding scores."),
    bbb: bool = typer.Option(False, "--bbb", help="Require BBB penetration prediction."),
) -> None:
    """Run safety and property screening on peptide candidates."""
    typer.echo(f"Screening not yet implemented. File: {candidates_file}")
