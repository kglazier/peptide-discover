"""Logging setup with Rich."""

import logging

from rich.console import Console
from rich.logging import RichHandler

console = Console()


def setup_logging(verbosity: int = 0) -> None:
    """Configure logging with Rich handler."""
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG

    logging.basicConfig(
        level=level,
        format="%(message)s",
        handlers=[RichHandler(console=console, rich_tracebacks=True)],
    )
