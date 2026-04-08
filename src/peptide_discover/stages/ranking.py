"""Stage 5: Multi-objective ranking and export."""

import json
import logging
from pathlib import Path

import pandas as pd

from peptide_discover.models.peptide import PeptideCandidate
from peptide_discover.models.scores import BindingResult, RankedCandidate, SafetyProfile

logger = logging.getLogger(__name__)


def rank_candidates(
    binding_results: list[BindingResult],
    safety_results: list[SafetyProfile],
    peptides: list[PeptideCandidate] | None = None,
    output_dir: Path = Path("results"),
) -> list[RankedCandidate]:
    """Rank candidates using direct scores from upstream tools.

    Design principle (from AMR project lesson): use simple, interpretable
    scores — not embedding-based meta-scoring. The neural networks already
    did the hard work inside PepMLM, Boltz-2, and PeptiVerse.

    Ranking is primarily by binding affinity (more negative = better).
    Safety scores are used as filters, not weights.
    """
    # Build lookup maps
    safety_map: dict[str, SafetyProfile] = {}
    for s in safety_results:
        safety_map[s.peptide_sequence] = s

    peptide_map: dict[str, PeptideCandidate] = {}
    if peptides:
        for p in peptides:
            peptide_map[p.sequence] = p

    # Filter out candidates flagged as toxic (if screening was run)
    filtered = []
    for br in binding_results:
        safety = safety_map.get(br.peptide_sequence)
        if safety and safety.toxicity_score is not None and safety.toxicity_score > 0.8:
            logger.info("Filtered %s: toxicity=%.2f", br.peptide_sequence, safety.toxicity_score)
            continue
        filtered.append(br)

    # Sort by affinity (more negative = better binding)
    filtered.sort(key=lambda r: r.affinity_score)

    # Build ranked candidates
    ranked = []
    for i, br in enumerate(filtered):
        peptide = peptide_map.get(
            br.peptide_sequence,
            PeptideCandidate(sequence=br.peptide_sequence),
        )
        safety = safety_map.get(
            br.peptide_sequence,
            SafetyProfile(peptide_sequence=br.peptide_sequence),
        )

        # Composite score: normalized affinity (primary signal)
        # Scale: -15 kcal/mol → 1.0, 0 kcal/mol → 0.0
        composite = min(1.0, max(0.0, abs(br.affinity_score) / 15.0))

        ranked.append(
            RankedCandidate(
                rank=i + 1,
                peptide=peptide,
                binding=br,
                safety=safety,
                composite_score=round(composite, 3),
            )
        )

    # Export
    output_dir.mkdir(parents=True, exist_ok=True)
    export_csv(ranked, output_dir / "ranked_candidates.csv")
    export_json(ranked, output_dir / "ranked_candidates.json")

    logger.info("Ranked %d candidates. Results in %s", len(ranked), output_dir)
    return ranked


def export_csv(candidates: list[RankedCandidate], output_path: Path) -> None:
    """Export ranked candidates to CSV."""
    rows = []
    for c in candidates:
        rows.append({
            "rank": c.rank,
            "sequence": c.peptide.sequence,
            "length": c.peptide.length,
            "generation_method": c.peptide.generation_method,
            "perplexity": c.peptide.perplexity,
            "affinity_kcal_mol": c.binding.affinity_score,
            "binding_confidence": c.binding.confidence,
            "toxicity": c.safety.toxicity_score,
            "bbb_penetration": c.safety.bbb_penetration,
            "solubility": c.safety.solubility,
            "half_life_hours": c.safety.half_life_hours,
            "hemolysis_risk": c.safety.hemolysis_risk,
            "composite_score": c.composite_score,
        })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    logger.info("Exported CSV: %s (%d rows)", output_path, len(rows))


def export_json(candidates: list[RankedCandidate], output_path: Path) -> None:
    """Export ranked candidates to JSON."""
    data = [c.model_dump() for c in candidates]
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2, default=str)
    logger.info("Exported JSON: %s", output_path)
