"""Scoring and ranking data models."""

from pydantic import BaseModel

from peptide_discover.models.peptide import PeptideCandidate


class BindingResult(BaseModel):
    """Binding prediction for a peptide-target pair."""

    peptide_sequence: str
    target_id: str
    affinity_score: float = 0.0
    confidence: float = 0.0
    complex_pdb_path: str | None = None


class SafetyProfile(BaseModel):
    """Safety/ADMET predictions for a peptide."""

    peptide_sequence: str
    toxicity_score: float | None = None
    bbb_penetration: bool | None = None
    solubility: float | None = None
    half_life_hours: float | None = None
    hemolysis_risk: float | None = None


class RankedCandidate(BaseModel):
    """A fully scored and ranked peptide candidate."""

    rank: int
    peptide: PeptideCandidate
    binding: BindingResult
    safety: SafetyProfile
    composite_score: float = 0.0
