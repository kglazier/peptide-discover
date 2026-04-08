"""Tests for safety screening stage."""

from peptide_discover.models.peptide import PeptideCandidate
from peptide_discover.stages.screening import screen_candidates, _compute_properties


def test_compute_properties():
    props = _compute_properties("ACDWK")
    assert props["length"] == 5
    assert props["mw"] > 0
    assert props["frac_cys"] == 0.2  # 1 C out of 5


def test_cysteine_rich_flagged_toxic():
    candidates = [PeptideCandidate(sequence="CCRCCRCCRC")]
    results = screen_candidates(candidates)
    assert results[0].toxicity_score > 0.5


def test_anionic_not_bbb():
    candidates = [PeptideCandidate(sequence="DDDDDDDDD")]
    results = screen_candidates(candidates, require_bbb=True)
    assert results[0].bbb_penetration is False


def test_bbb_not_computed_when_not_required():
    candidates = [PeptideCandidate(sequence="ACDWK")]
    results = screen_candidates(candidates, require_bbb=False)
    assert results[0].bbb_penetration is None


def test_hydrophobic_low_solubility():
    candidates = [PeptideCandidate(sequence="LLLLLLLLLL")]
    results = screen_candidates(candidates)
    assert results[0].solubility < 0.2


def test_charged_high_solubility():
    candidates = [PeptideCandidate(sequence="RKRKRKRKRK")]
    results = screen_candidates(candidates)
    assert results[0].solubility > 0.8


def test_batch_screening():
    candidates = [
        PeptideCandidate(sequence="ACDWK"),
        PeptideCandidate(sequence="RLHFW"),
        PeptideCandidate(sequence="GGGGG"),
    ]
    results = screen_candidates(candidates)
    assert len(results) == 3
    assert all(r.peptide_sequence for r in results)
