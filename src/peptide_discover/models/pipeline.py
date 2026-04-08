"""Full pipeline result model."""

from datetime import datetime

from pydantic import BaseModel

from peptide_discover.models.scores import RankedCandidate
from peptide_discover.models.target import TargetProtein


class PipelineResult(BaseModel):
    """Complete output from a discovery pipeline run."""

    target: TargetProtein
    track: str = ""
    candidates: list[RankedCandidate] = []
    total_generated: int = 0
    total_after_binding_filter: int = 0
    started_at: datetime = datetime.now()
    completed_at: datetime | None = None
