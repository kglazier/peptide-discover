"""BPC-157 docking against gut-relevant targets.

Docks BPC-157 (GEPPPGKPADDAGLV) against targets relevant to its
claimed mechanisms: angiogenesis, gut repair, wound healing, NO signaling.

Usage:
    python experiments/bpc157_docking.py
"""

import json
import logging
from pathlib import Path

from peptide_discover.models.peptide import PeptideCandidate
from peptide_discover.stages.binding import predict_binding
from peptide_discover.stages.screening import screen_candidates
from peptide_discover.stages.target_prep import resolve_target

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

BPC157 = "GEPPPGKPADDAGLV"

TARGETS = {
    "VEGFR2": ("P35968", "Angiogenesis — new blood vessel formation"),
    "EGFR": ("P00533", "Epithelial growth factor receptor — gut lining repair"),
    "PDGFRB": ("P09619", "Wound healing / tissue repair"),
    "NOS3": ("P29474", "Endothelial nitric oxide synthase — NO signaling"),
}

RESULTS_DIR = Path("results/bpc157")


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    candidate = PeptideCandidate(sequence=BPC157, generation_method="known")

    # Safety screening
    logger.info("Safety screening BPC-157...")
    safety = screen_candidates([candidate], require_bbb=False)
    s = safety[0]
    logger.info("  Toxicity: %s, Solubility: %s, Hemolysis: %s",
                s.toxicity_score, s.solubility, s.hemolysis_risk)

    # Dock against each target
    results = []
    for name, (uniprot, role) in TARGETS.items():
        logger.info("Target: %s (%s) — %s", name, uniprot, role)
        try:
            target = resolve_target(uniprot, output_dir=RESULTS_DIR / "targets")
            binding = predict_binding(target, [candidate], top_k=1, exhaustiveness=4)
            score = binding[0].affinity_score if binding else 0.0
            logger.info("  %s: %.2f kcal/mol", name, score)
        except Exception as e:
            score = 0.0
            logger.error("  %s: error — %s", name, str(e)[:80])

        results.append({
            "target": name,
            "uniprot": uniprot,
            "role": role,
            "affinity_kcal_mol": score,
        })

    # Summary
    print(f"\n{'Target':<10} {'Role':<50} {'Affinity':>10}")
    print("-" * 72)
    for r in results:
        print(f"{r['target']:<10} {r['role']:<50} {r['affinity_kcal_mol']:>8.2f}")

    # Save
    output = {
        "peptide": "BPC-157",
        "sequence": BPC157,
        "safety": {
            "toxicity": s.toxicity_score,
            "solubility": s.solubility,
            "hemolysis_risk": s.hemolysis_risk,
        },
        "docking_results": results,
    }
    with open(RESULTS_DIR / "bpc157_results.json", "w") as f:
        json.dump(output, f, indent=2)
    logger.info("Saved to %s", RESULTS_DIR / "bpc157_results.json")


if __name__ == "__main__":
    main()
