"""Test Pydantic data models."""

from peptide_discover.models.peptide import PeptideCandidate
from peptide_discover.models.scores import BindingResult, RankedCandidate, SafetyProfile
from peptide_discover.models.target import BindingSite, TargetProtein


def test_peptide_auto_length():
    p = PeptideCandidate(sequence="ACDWKTHLYG")
    assert p.length == 10


def test_peptide_explicit_length():
    p = PeptideCandidate(sequence="ACDWKTHLYG", length=10)
    assert p.length == 10


def test_target_protein():
    t = TargetProtein(identifier="P04629", sequence="MLRG")
    assert t.identifier == "P04629"
    assert t.structure_path is None


def test_binding_site():
    bs = BindingSite(chain="A", residue_indices=[10, 20, 30])
    assert len(bs.residue_indices) == 3


def test_binding_result():
    br = BindingResult(
        peptide_sequence="ACDWK",
        target_id="P04629",
        affinity_score=-8.5,
        confidence=0.92,
    )
    assert br.affinity_score == -8.5


def test_safety_profile():
    sp = SafetyProfile(
        peptide_sequence="ACDWK",
        toxicity_score=0.1,
        bbb_penetration=True,
    )
    assert sp.bbb_penetration is True


def test_ranked_candidate():
    rc = RankedCandidate(
        rank=1,
        peptide=PeptideCandidate(sequence="ACDWK"),
        binding=BindingResult(peptide_sequence="ACDWK", target_id="P04629"),
        safety=SafetyProfile(peptide_sequence="ACDWK"),
        composite_score=0.85,
    )
    assert rc.rank == 1
    assert rc.composite_score == 0.85
