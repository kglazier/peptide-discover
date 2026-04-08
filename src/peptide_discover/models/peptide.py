"""Peptide candidate data models."""

from pydantic import BaseModel


class PeptideCandidate(BaseModel):
    """A generated peptide candidate."""

    sequence: str
    length: int = 0
    generation_method: str = ""  # pepmlm, rfdiffusion, etc.
    generation_rank: int = 0     # rank from the generation model
    perplexity: float | None = None  # pseudo-perplexity score (lower = better)

    def model_post_init(self, __context) -> None:
        if self.length == 0:
            self.length = len(self.sequence)
