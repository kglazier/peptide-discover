"""Full pipeline orchestration command."""

from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

app = typer.Typer()
console = Console()


@app.command("pipeline")
def run_pipeline(
    target: str = typer.Argument(help="Target protein: UniProt ID, PDB ID, or file path."),
    track: Optional[str] = typer.Option(  # noqa: UP007
        None, "--track", "-t", help="Research track: cognitive, muscle, or amp."
    ),
    candidates: int = typer.Option(100, "--candidates", "-n"),
    peptide_length: int = typer.Option(15, "--length", "-l", help="Peptide length (3-50)."),
    top_k_gen: int = typer.Option(3, "--top-k-gen", help="Top-k sampling for generation."),
    method: str = typer.Option("pepmlm", "--method", "-m", help="Generation method."),
    output_dir: Path = typer.Option(Path("results"), "--output-dir", "-o"),
    skip_screening: bool = typer.Option(False, "--skip-screening"),
    skip_docking: bool = typer.Option(False, "--skip-docking"),
    top_k: int = typer.Option(50, "--top-k", "-k", help="Keep top K after docking."),
    engine: str = typer.Option("vina", "--engine", help="Docking engine: vina (fast, 8-mers) or adcp (handles 15+ mers, WSL)."),
    exhaustiveness: int = typer.Option(8, "--exhaustiveness", "-e", help="Vina exhaustiveness."),
    adcp_runs: int = typer.Option(10, "--adcp-runs", help="ADCP MC runs per peptide."),
    adcp_steps: int = typer.Option(2500000, "--adcp-steps", help="ADCP MC steps per run."),
) -> None:
    """Run the full peptide discovery pipeline."""
    from peptide_discover.config.tracks import get_track
    from peptide_discover.stages.target_prep import resolve_target
    from peptide_discover.stages.generation import generate_pepmlm
    from peptide_discover.stages.binding import predict_binding
    from peptide_discover.stages.ranking import rank_candidates

    # Resolve track if provided
    research_track = get_track(track) if track else None

    console.print()
    console.print("[bold]peptide-discover[/bold] pipeline")
    console.print(f"  Target:     {target}")
    if research_track:
        console.print(f"  Track:      {research_track.name} — {research_track.description}")
    console.print(f"  Candidates: {candidates}")
    console.print(f"  Length:     {peptide_length}aa")
    console.print(f"  Method:     {method}")
    console.print()

    # Stage 1: Target preparation
    console.print("[bold cyan]Stage 1:[/bold cyan] Preparing target...")
    target_protein = resolve_target(target, output_dir=output_dir / "targets")
    console.print(f"  {target_protein.identifier} — {len(target_protein.sequence)}aa")
    if target_protein.structure_path:
        console.print(f"  Structure: {target_protein.structure_path}")
    console.print()

    # Stage 2: Peptide generation
    console.print(f"[bold cyan]Stage 2:[/bold cyan] Generating {candidates} candidates...")
    peptides = generate_pepmlm(
        target_protein, n=candidates, peptide_length=peptide_length, top_k=top_k_gen,
    )
    console.print(f"  {len(peptides)} unique candidates generated.")
    console.print()

    # Stage 3: Binding prediction (optional)
    if not skip_docking:
        console.print(f"[bold cyan]Stage 3:[/bold cyan] Docking peptides with {engine} (this may take a while)...")
        binding_results = predict_binding(
            target_protein, peptides, top_k=top_k,
            engine=engine, exhaustiveness=exhaustiveness,
            adcp_runs=adcp_runs, adcp_steps=adcp_steps,
        )
        console.print(f"  {len(binding_results)} candidates scored.")
    else:
        console.print("[dim]Stage 3: Skipped (--skip-docking).[/dim]")
        # Create placeholder binding results sorted by perplexity
        from peptide_discover.models.scores import BindingResult
        binding_results = [
            BindingResult(
                peptide_sequence=p.sequence,
                target_id=target_protein.identifier,
                affinity_score=0.0,
                confidence=0.0,
            )
            for p in peptides
        ]
    console.print()

    # Stage 4: Safety screening (optional, not yet implemented)
    safety_results = []
    if not skip_screening:
        try:
            from peptide_discover.stages.screening import screen_candidates
            require_bbb = research_track.require_bbb if research_track else False
            console.print("[bold cyan]Stage 4:[/bold cyan] Running safety screening...")
            safety_results = screen_candidates(peptides, require_bbb=require_bbb)
        except NotImplementedError:
            console.print("[dim]Stage 4: Safety screening not yet implemented, skipping.[/dim]")
    else:
        console.print("[dim]Stage 4: Skipped (--skip-screening).[/dim]")
    console.print()

    # Stage 5: Ranking & export
    console.print("[bold cyan]Stage 5:[/bold cyan] Ranking and exporting...")
    ranked = rank_candidates(
        binding_results, safety_results, peptides=peptides, output_dir=output_dir,
    )

    # Display results table
    table = Table(title=f"Top {min(20, len(ranked))} candidates")
    table.add_column("Rank", justify="right", style="cyan")
    table.add_column("Sequence", style="bold")
    table.add_column("Affinity", justify="right", style="green")
    table.add_column("Perplexity", justify="right")
    table.add_column("Score", justify="right", style="yellow")

    for c in ranked[:20]:
        aff = f"{c.binding.affinity_score:.1f}" if c.binding.affinity_score != 0 else "—"
        ppl = f"{c.peptide.perplexity:.1f}" if c.peptide.perplexity else "—"
        table.add_row(
            str(c.rank),
            c.peptide.sequence,
            aff,
            ppl,
            f"{c.composite_score:.3f}",
        )
    console.print(table)

    console.print()
    console.print(f"[bold green]Done.[/bold green] {len(ranked)} candidates ranked.")
    console.print(f"  CSV:  {output_dir / 'ranked_candidates.csv'}")
    console.print(f"  JSON: {output_dir / 'ranked_candidates.json'}")
