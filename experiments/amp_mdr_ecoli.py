"""AMP design for multi-drug resistant E. coli.

Cross-project experiment linking peptide-discover and AMRCast:
1. Analyze AMRCast resistance data to identify MDR patterns
2. Design antimicrobial peptides targeting E. coli membrane proteins
3. Dock and screen candidates
4. Report: novel AMP candidates for MDR E. coli strains

Usage:
    python experiments/amp_mdr_ecoli.py [--skip-docking] [--candidates 20]

Requires:
    - peptide-discover (this repo)
    - AMRCast data at ../biology/amrcast/data/narms/
"""

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
AMRCAST_DATA = PROJECT_ROOT.parent / "biology" / "amrcast" / "data" / "narms"
RESULTS_DIR = PROJECT_ROOT / "results" / "amp_mdr"

# E. coli outer membrane protein targets for AMPs
AMP_TARGETS = {
    "OmpA": {
        "uniprot": "P0A910",
        "role": "Outer membrane protein A — major surface protein, structural integrity",
    },
    "BamA": {
        "uniprot": "P0A940",
        "role": "Outer membrane assembly factor — essential for membrane biogenesis",
    },
    "LptD": {
        "uniprot": "P31554",
        "role": "LPS transport — essential for outer membrane LPS assembly",
    },
}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def analyze_resistance(antibiogram_path: Path) -> dict:
    """Analyze AMRCast resistance data to understand the MDR landscape."""
    logger.info("Loading AMRCast antibiogram data...")
    df = pd.read_csv(antibiogram_path)

    total_genomes = df.biosample_acc.nunique()
    resistant = df[df.resistance_phenotype == "resistant"]

    # MDR = resistant to 3+ drug classes
    mdr = resistant.groupby("biosample_acc").antibiotic.nunique()
    mdr_counts = {
        "total_genomes": total_genomes,
        "any_resistance": resistant.biosample_acc.nunique(),
        "mdr_3plus": int((mdr >= 3).sum()),
        "mdr_5plus": int((mdr >= 5).sum()),
        "mdr_10plus": int((mdr >= 10).sum()),
    }

    # Which antibiotics have the most resistance?
    resistance_rates = (
        resistant.groupby("antibiotic")
        .biosample_acc.nunique()
        .sort_values(ascending=False)
    )
    top_resistant = resistance_rates.head(10)

    logger.info("Resistance landscape:")
    logger.info("  Total genomes: %d", total_genomes)
    logger.info("  Any resistance: %d (%.1f%%)", mdr_counts["any_resistance"],
                100 * mdr_counts["any_resistance"] / total_genomes)
    logger.info("  MDR (5+ drugs): %d (%.1f%%)", mdr_counts["mdr_5plus"],
                100 * mdr_counts["mdr_5plus"] / total_genomes)
    logger.info("  Top resistant antibiotics:")
    for ab, count in top_resistant.items():
        logger.info("    %s: %d strains (%.1f%%)", ab, count, 100 * count / total_genomes)

    return {
        "counts": mdr_counts,
        "top_resistant_antibiotics": {
            ab: {"strains": int(count), "pct": round(100 * count / total_genomes, 1)}
            for ab, count in top_resistant.items()
        },
    }


def design_amps(n_candidates: int, peptide_length: int, skip_docking: bool, exhaustiveness: int = 4) -> dict:
    """Design AMPs against E. coli membrane targets."""
    from peptide_discover.stages.target_prep import resolve_target
    from peptide_discover.stages.generation import generate_pepmlm
    from peptide_discover.stages.binding import predict_binding
    from peptide_discover.stages.screening import screen_candidates
    from peptide_discover.stages.ranking import rank_candidates

    all_results = {}

    for name, info in AMP_TARGETS.items():
        logger.info("=" * 60)
        logger.info("Target: %s — %s", name, info["role"])
        logger.info("UniProt: %s", info["uniprot"])

        # Stage 1: Target prep
        target = resolve_target(info["uniprot"], output_dir=RESULTS_DIR / "targets")
        logger.info("  Sequence: %d aa", len(target.sequence))

        # Stage 2: Generate candidates
        logger.info("  Generating %d candidates...", n_candidates)
        peptides = generate_pepmlm(
            target, n=n_candidates, peptide_length=peptide_length, top_k=5,
        )
        logger.info("  Generated %d unique peptides.", len(peptides))

        # Stage 3: Binding prediction
        if not skip_docking:
            logger.info("  Docking (this will take a while)...")
            binding_results = predict_binding(
                target, peptides, top_k=n_candidates, exhaustiveness=exhaustiveness,
            )
        else:
            logger.info("  Docking skipped.")
            from peptide_discover.models.scores import BindingResult
            binding_results = [
                BindingResult(
                    peptide_sequence=p.sequence,
                    target_id=target.identifier,
                )
                for p in peptides
            ]

        # Stage 4: Safety screening (AMPs need BBB=False, but do need toxicity)
        logger.info("  Screening for safety...")
        safety_results = screen_candidates(peptides, require_bbb=False)

        # Stage 5: Rank
        ranked = rank_candidates(
            binding_results, safety_results, peptides=peptides,
            output_dir=RESULTS_DIR / name,
        )

        all_results[name] = {
            "target_uniprot": info["uniprot"],
            "target_role": info["role"],
            "candidates_generated": len(peptides),
            "candidates_after_ranking": len(ranked),
            "top_5": [
                {
                    "rank": c.rank,
                    "sequence": c.peptide.sequence,
                    "perplexity": c.peptide.perplexity,
                    "affinity_kcal_mol": c.binding.affinity_score,
                    "toxicity": c.safety.toxicity_score,
                    "hemolysis_risk": c.safety.hemolysis_risk,
                    "solubility": c.safety.solubility,
                    "composite_score": c.composite_score,
                }
                for c in ranked[:5]
            ],
        }

        logger.info("  Top candidates for %s:", name)
        for c in ranked[:5]:
            aff = f"{c.binding.affinity_score:.1f}" if c.binding.affinity_score else "—"
            logger.info(
                "    #%d %s  affinity=%s  tox=%.2f  sol=%.2f",
                c.rank, c.peptide.sequence, aff,
                c.safety.toxicity_score or 0, c.safety.solubility or 0,
            )

    return all_results


def generate_report(resistance_data: dict, amp_data: dict, output_path: Path) -> None:
    """Generate a combined report."""
    report = {
        "title": "AMP Design for Multi-Drug Resistant E. coli",
        "generated": datetime.now().isoformat(),
        "motivation": (
            "With {mdr5} E. coli strains ({mdr5_pct:.1f}%) resistant to 5+ antibiotics, "
            "antimicrobial peptides offer an alternative treatment mechanism that "
            "bypasses conventional resistance pathways."
        ).format(
            mdr5=resistance_data["counts"]["mdr_5plus"],
            mdr5_pct=100 * resistance_data["counts"]["mdr_5plus"]
            / resistance_data["counts"]["total_genomes"],
        ),
        "resistance_landscape": resistance_data,
        "amp_candidates": amp_data,
        "methodology": {
            "target_selection": "E. coli outer membrane proteins essential for viability",
            "generation": "PepMLM (ESM-2 650M fine-tuned for peptide binder design)",
            "binding": "AutoDock Vina molecular docking",
            "screening": "Local physicochemical property-based safety scoring",
            "data_source": "AMRCast NARMS antibiogram database",
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(report, f, indent=2, default=str)
    logger.info("Report saved: %s", output_path)


def main():
    parser = argparse.ArgumentParser(description="AMP design for MDR E. coli")
    parser.add_argument("--candidates", "-n", type=int, default=30,
                        help="Candidates per target (default: 30)")
    parser.add_argument("--length", "-l", type=int, default=8,
                        help="Peptide length (default: 8)")
    parser.add_argument("--exhaustiveness", "-e", type=int, default=4,
                        help="Vina exhaustiveness (default: 4)")
    parser.add_argument("--skip-docking", action="store_true",
                        help="Skip docking (fast mode for testing)")
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Check AMRCast data is available
    antibiogram = AMRCAST_DATA / "antibiogram_mic.csv"
    if not antibiogram.exists():
        logger.error("AMRCast data not found at %s", antibiogram)
        logger.error("Expected: projects/biology/amrcast/data/narms/antibiogram_mic.csv")
        sys.exit(1)

    # Step 1: Analyze resistance landscape
    logger.info("Step 1: Analyzing resistance landscape from AMRCast data...")
    resistance_data = analyze_resistance(antibiogram)

    # Step 2: Design AMPs against E. coli membrane targets
    logger.info("\nStep 2: Designing AMPs against E. coli membrane targets...")
    amp_data = design_amps(
        n_candidates=args.candidates,
        peptide_length=args.length,
        skip_docking=args.skip_docking,
        exhaustiveness=args.exhaustiveness,
    )

    # Step 3: Generate combined report
    logger.info("\nStep 3: Generating report...")
    generate_report(
        resistance_data, amp_data,
        RESULTS_DIR / "amp_mdr_ecoli_report.json",
    )

    logger.info("\nDone. Results in %s", RESULTS_DIR)


if __name__ == "__main__":
    main()
