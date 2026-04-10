"""Baseline experiment: dock known peptides against their real targets.

Sanity check for the pipeline — if our docking produces reasonable scores
for peptides that are known to bind their targets, we can trust the scores
we get for novel candidates.

Usage:
    python experiments/baseline_known_peptides.py
"""

import logging
from pathlib import Path

from peptide_discover.models.peptide import PeptideCandidate
from peptide_discover.stages.binding import predict_binding
from peptide_discover.stages.target_prep import resolve_target

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# Known peptides with their real targets
# Each entry: (peptide_name, sequence, target_uniprot, target_name, notes)
KNOWN_PEPTIDES = [
    # Cognitive / neuroplasticity
    (
        "Semax",
        "MEHFPGP",
        "P01189",  # POMC (ACTH) — Semax is a fragment analog of ACTH 4-10
        "POMC",
        "Nootropic, neuroprotection",
    ),
    (
        "Selank",
        "TKPRPGP",
        "P01019",  # AGT (tuftsin-derived) — Selank binds GABA receptors
        "AGT",
        "Anxiolytic nootropic",
    ),
    # Muscle / GH secretagogues
    (
        "Ipamorelin",
        "AIWFK",  # 5-mer simplified, actual has D-amino acids
        "Q92847",  # GHSR
        "GHSR",
        "GH secretagogue receptor agonist",
    ),
    # Antimicrobial peptides (natural AMPs)
    (
        "Magainin-2-short",
        "GIGKFLHS",  # 8-mer fragment of full 23-mer magainin
        "P0A910",  # OmpA
        "OmpA",
        "Natural AMP from Xenopus skin",
    ),
    (
        "LL-37-short",
        "LLGDFFRK",  # 8-mer fragment of full LL-37
        "P0A910",  # OmpA
        "OmpA",
        "Human cathelicidin fragment",
    ),
    (
        "Cecropin-A-short",
        "KWKLFKKI",  # 8-mer fragment of full cecropin
        "P0A910",  # OmpA
        "OmpA",
        "Insect AMP",
    ),
    (
        "Magainin-2-bamA",
        "GIGKFLHS",
        "P0A940",  # BamA
        "BamA",
        "Natural AMP vs BamA",
    ),
    (
        "Magainin-2-lptD",
        "GIGKFLHS",
        "P31554",  # LptD
        "LptD",
        "Natural AMP vs LptD",
    ),
]


def main():
    results = []

    # Group by target to avoid re-docking the same receptor
    targets_seen = {}

    for name, sequence, uniprot, target_name, notes in KNOWN_PEPTIDES:
        logger.info("=" * 60)
        logger.info("%s (%s) vs %s (%s)", name, sequence, target_name, uniprot)
        logger.info("  %s", notes)

        # Resolve target (cached after first fetch)
        if uniprot not in targets_seen:
            targets_seen[uniprot] = resolve_target(uniprot)
        target = targets_seen[uniprot]

        # Dock the single peptide
        candidate = PeptideCandidate(
            sequence=sequence,
            generation_method="known",
            generation_rank=1,
        )

        try:
            binding_results = predict_binding(
                target, [candidate], top_k=1, exhaustiveness=4,
            )
            if binding_results:
                score = binding_results[0].affinity_score
                logger.info("  Affinity: %.2f kcal/mol", score)
                results.append({
                    "peptide": name,
                    "sequence": sequence,
                    "target": target_name,
                    "affinity_kcal_mol": score,
                    "notes": notes,
                })
            else:
                logger.warning("  Failed to dock")
        except Exception as e:
            logger.error("  Error: %s", str(e)[:100])

    # Print summary
    logger.info("\n" + "=" * 60)
    logger.info("BASELINE SUMMARY")
    logger.info("=" * 60)
    print(f"\n{'Peptide':<20} {'Sequence':<12} {'Target':<8} {'Affinity':>10}")
    print("-" * 54)
    for r in results:
        print(f"{r['peptide']:<20} {r['sequence']:<12} {r['target']:<8} {r['affinity_kcal_mol']:>8.2f}")

    # Save results
    import json
    output_dir = Path("results/baseline")
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "known_peptides.json", "w") as f:
        json.dump(results, f, indent=2)
    logger.info("\nSaved to %s", output_dir / "known_peptides.json")


if __name__ == "__main__":
    main()
