"""Stage 4: Property and safety screening.

Local physicochemical property-based screening. Computes well-established
peptide properties and flags candidates based on known thresholds from
the literature. Not as accurate as a trained ML model, but transparent,
fast, and doesn't depend on external services.

References:
- Gupta et al., ToxinPred 3.0 (2024) — AAC/DPC features for toxicity
- Kumar et al., B3Pred (2021) — physicochemical features for BBB penetration
- Kyte & Doolittle (1982) — hydrophobicity scale
"""

import logging

from peptide_discover.models.peptide import PeptideCandidate
from peptide_discover.models.scores import SafetyProfile

logger = logging.getLogger(__name__)

# Kyte-Doolittle hydrophobicity scale
HYDROPHOBICITY = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}

# Amino acid molecular weights (Da)
AA_MW = {
    "A": 89.1, "C": 121.2, "D": 133.1, "E": 147.1, "F": 165.2,
    "G": 75.0, "H": 155.2, "I": 131.2, "K": 146.2, "L": 131.2,
    "M": 149.2, "N": 132.1, "P": 115.1, "Q": 146.2, "R": 174.2,
    "S": 105.1, "T": 119.1, "V": 117.1, "W": 204.2, "Y": 181.2,
}

# Amino acid charge at pH 7.4
AA_CHARGE = {
    "D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1,
}

# Hydrophobic residues
HYDROPHOBIC_AAS = set("AILMFVWP")

# Cationic residues
CATIONIC_AAS = set("KRH")


def screen_candidates(
    candidates: list[PeptideCandidate],
    require_bbb: bool = False,
) -> list[SafetyProfile]:
    """Run safety screening using local physicochemical property analysis.

    Computes:
    - Toxicity score: based on cysteine content, charge, and toxic motif patterns
    - BBB penetration: based on charge, hydrophobicity, amphipathicity, and MW
    - Solubility estimate: based on charge and hydrophobicity balance
    - Hemolysis risk: based on hydrophobicity and cationic character
    """
    profiles = []
    for candidate in candidates:
        seq = candidate.sequence
        props = _compute_properties(seq)

        toxicity = _score_toxicity(seq, props)
        bbb = _predict_bbb(seq, props) if require_bbb else None
        solubility = _score_solubility(props)
        hemolysis = _score_hemolysis(seq, props)

        profiles.append(
            SafetyProfile(
                peptide_sequence=seq,
                toxicity_score=toxicity,
                bbb_penetration=bbb,
                solubility=solubility,
                hemolysis_risk=hemolysis,
            )
        )

    toxic_count = sum(1 for p in profiles if p.toxicity_score and p.toxicity_score > 0.5)
    bbb_count = sum(1 for p in profiles if p.bbb_penetration is True)
    logger.info(
        "Screening complete. %d/%d flagged toxic. %d/%d predicted BBB-penetrating.",
        toxic_count, len(profiles),
        bbb_count, len(profiles),
    )
    return profiles


def _compute_properties(seq: str) -> dict:
    """Compute physicochemical properties for a peptide sequence."""
    length = len(seq)
    if length == 0:
        return {}

    # Molecular weight
    mw = sum(AA_MW.get(aa, 110.0) for aa in seq) - 18.02 * (length - 1)

    # Net charge at pH 7.4
    charge = sum(AA_CHARGE.get(aa, 0.0) for aa in seq)

    # Mean hydrophobicity
    hydro_values = [HYDROPHOBICITY.get(aa, 0.0) for aa in seq]
    mean_hydro = sum(hydro_values) / length

    # Fraction hydrophobic
    frac_hydrophobic = sum(1 for aa in seq if aa in HYDROPHOBIC_AAS) / length

    # Fraction cationic
    frac_cationic = sum(1 for aa in seq if aa in CATIONIC_AAS) / length

    # Cysteine fraction
    frac_cys = seq.count("C") / length

    # Amphipathicity estimate: variance in hydrophobicity
    # High variance = amphipathic (mix of hydrophobic and hydrophilic)
    if length > 1:
        amphipathicity = sum((h - mean_hydro) ** 2 for h in hydro_values) / length
    else:
        amphipathicity = 0.0

    return {
        "length": length,
        "mw": mw,
        "charge": charge,
        "mean_hydro": mean_hydro,
        "frac_hydrophobic": frac_hydrophobic,
        "frac_cationic": frac_cationic,
        "frac_cys": frac_cys,
        "amphipathicity": amphipathicity,
    }


def _score_toxicity(seq: str, props: dict) -> float:
    """Score toxicity risk (0-1, higher = more likely toxic).

    Based on features from ToxinPred literature:
    - Animal toxins are typically cysteine-rich (disulfide bonds for stability)
    - Many toxins are highly cationic
    - Specific dipeptide patterns associated with toxicity
    """
    score = 0.0

    # Cysteine content: toxins average ~12% Cys, non-toxins ~1.5%
    frac_cys = props.get("frac_cys", 0)
    if frac_cys > 0.10:
        score += 0.35
    elif frac_cys > 0.05:
        score += 0.15

    # High positive charge: many toxins are strongly cationic
    charge = props.get("charge", 0)
    charge_per_residue = charge / max(props.get("length", 1), 1)
    if charge_per_residue > 0.3:
        score += 0.20
    elif charge_per_residue > 0.2:
        score += 0.10

    # Very hydrophobic peptides can be membrane-disrupting (toxic)
    if props.get("mean_hydro", 0) > 2.0:
        score += 0.15

    # Known toxic motifs: CxC patterns (disulfide frameworks)
    import re
    cxc_count = len(re.findall(r"C.{1,6}C", seq))
    if cxc_count >= 2:
        score += 0.25
    elif cxc_count == 1:
        score += 0.10

    # Very short peptides are generally less toxic
    if props.get("length", 0) < 6:
        score *= 0.5

    return min(1.0, round(score, 3))


def _predict_bbb(seq: str, props: dict) -> bool:
    """Predict blood-brain barrier penetration.

    BBB-penetrating peptides tend to be:
    - Short (5-30 residues, ideal 5-20)
    - Cationic (net positive charge)
    - Amphipathic (mix of hydrophobic and hydrophilic)
    - Low molecular weight (<4000 Da)
    - Moderate hydrophobicity

    Based on features used in B3Pred (Kumar et al., 2021).
    """
    length = props.get("length", 0)
    mw = props.get("mw", 0)
    charge = props.get("charge", 0)
    mean_hydro = props.get("mean_hydro", 0)
    amphipathicity = props.get("amphipathicity", 0)
    frac_cationic = props.get("frac_cationic", 0)

    bbb_score = 0.0

    # Length: 5-20 is ideal
    if 5 <= length <= 20:
        bbb_score += 1.0
    elif 20 < length <= 30:
        bbb_score += 0.5

    # MW: <4000 Da favored
    if mw < 2500:
        bbb_score += 1.0
    elif mw < 4000:
        bbb_score += 0.5

    # Cationic: positive charge helps transcytosis
    if charge >= 2:
        bbb_score += 1.0
    elif charge >= 1:
        bbb_score += 0.5

    # Moderate hydrophobicity (not too hydrophilic, not too hydrophobic)
    if -0.5 <= mean_hydro <= 1.5:
        bbb_score += 0.5
    elif mean_hydro > 1.5:
        bbb_score += 0.3

    # Amphipathicity: higher is better for membrane interaction
    if amphipathicity > 5.0:
        bbb_score += 1.0
    elif amphipathicity > 2.0:
        bbb_score += 0.5

    # Threshold: 3.0 out of 4.5 max
    return bbb_score >= 3.0


def _score_solubility(props: dict) -> float:
    """Estimate aqueous solubility (0-1, higher = more soluble).

    Highly hydrophobic peptides tend to aggregate.
    Charged peptides tend to be more soluble.
    """
    mean_hydro = props.get("mean_hydro", 0)
    charge = abs(props.get("charge", 0))
    frac_hydrophobic = props.get("frac_hydrophobic", 0)

    # Start with inverse hydrophobicity
    # Scale: mean_hydro ranges from ~-4 to ~+4
    solubility = 0.5 - (mean_hydro / 8.0)  # -4 → 1.0, +4 → 0.0

    # Charge helps solubility
    solubility += min(0.3, charge * 0.1)

    # Very high hydrophobic fraction is bad
    if frac_hydrophobic > 0.7:
        solubility -= 0.2

    return max(0.0, min(1.0, round(solubility, 3)))


def _score_hemolysis(seq: str, props: dict) -> float:
    """Score hemolysis risk (0-1, higher = more hemolytic).

    Hemolytic peptides tend to be:
    - Amphipathic
    - Cationic
    - Hydrophobic
    """
    score = 0.0

    # Amphipathic + cationic = classic hemolytic profile
    if props.get("amphipathicity", 0) > 5.0 and props.get("frac_cationic", 0) > 0.2:
        score += 0.3

    # High hydrophobicity
    if props.get("mean_hydro", 0) > 1.0:
        score += 0.2

    # Leucine-rich (common in hemolytic peptides)
    leu_frac = seq.count("L") / max(len(seq), 1)
    if leu_frac > 0.3:
        score += 0.2

    # Strong positive charge
    if props.get("charge", 0) > 3:
        score += 0.15

    return min(1.0, round(score, 3))
