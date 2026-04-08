"""Root CLI application."""

from typing import Optional

import typer

from peptide_discover import __version__

app = typer.Typer(
    name="peptide-discover",
    help="AI peptide discovery pipeline.",
    no_args_is_help=True,
)


def version_callback(value: bool) -> None:
    if value:
        typer.echo(f"peptide-discover {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Optional[bool] = typer.Option(  # noqa: UP007
        None,
        "--version",
        "-V",
        help="Show version and exit.",
        callback=version_callback,
        is_eager=True,
    ),
) -> None:
    """AI peptide discovery pipeline."""


# Register subcommands
from peptide_discover.cli.target_cmd import app as target_app  # noqa: E402
from peptide_discover.cli.generate_cmd import app as generate_app  # noqa: E402
from peptide_discover.cli.binding_cmd import app as binding_app  # noqa: E402
from peptide_discover.cli.screen_cmd import app as screen_app  # noqa: E402
from peptide_discover.cli.rank_cmd import app as rank_app  # noqa: E402
from peptide_discover.cli.run_cmd import app as run_app  # noqa: E402

app.add_typer(target_app, name="target", help="Fetch and prepare target protein.")
app.add_typer(generate_app, name="generate", help="Generate peptide candidates.")
app.add_typer(binding_app, name="bind", help="Predict binding affinity.")
app.add_typer(screen_app, name="screen", help="Run property/safety screening.")
app.add_typer(rank_app, name="rank", help="Rank and export results.")
app.add_typer(run_app, name="run", help="Run full discovery pipeline.")
