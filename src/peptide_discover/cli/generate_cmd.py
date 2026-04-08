"""CLI commands for peptide generation."""

import json
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

app = typer.Typer()
console = Console()


@app.command("pepmlm")
def generate_pepmlm(
    target: str = typer.Argument(help="Target protein identifier or path."),
    candidates: int = typer.Option(100, "--candidates", "-n"),
    peptide_length: int = typer.Option(15, "--length", "-l", help="Peptide length (3-50)."),
    top_k: int = typer.Option(3, "--top-k", "-k", help="Top-k sampling (higher = more diverse)."),
    output_dir: Path = typer.Option(Path("results/generated"), "--output-dir", "-o"),
) -> None:
    """Generate peptide candidates using PepMLM."""
    from peptide_discover.stages.generation import generate_pepmlm as _generate
    from peptide_discover.stages.target_prep import resolve_target

    target_protein = resolve_target(target)
    console.print(f"Target: {target_protein.identifier} ({len(target_protein.sequence)}aa)")

    peptides = _generate(target_protein, n=candidates, peptide_length=peptide_length, top_k=top_k)

    # Display results
    table = Table(title=f"Generated {len(peptides)} peptide candidates")
    table.add_column("Rank", justify="right", style="cyan")
    table.add_column("Sequence", style="bold")
    table.add_column("Length", justify="right")
    table.add_column("Perplexity", justify="right", style="green")

    for p in peptides:
        table.add_row(
            str(p.generation_rank),
            p.sequence,
            str(p.length),
            f"{p.perplexity:.2f}" if p.perplexity else "—",
        )
    console.print(table)

    # Save to file
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{target_protein.identifier}_candidates.json"
    with open(output_file, "w") as f:
        json.dump([p.model_dump() for p in peptides], f, indent=2)
    console.print(f"Saved to {output_file}")
