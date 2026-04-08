"""Shared test fixtures."""

from pathlib import Path

import pytest

from peptide_discover.config.settings import Settings
from peptide_discover.models.peptide import PeptideCandidate
from peptide_discover.models.target import TargetProtein


@pytest.fixture
def test_data_dir() -> Path:
    return Path(__file__).parent / "data"


@pytest.fixture
def settings(tmp_path: Path) -> Settings:
    return Settings(
        data_dir=tmp_path / "data",
        results_dir=tmp_path / "results",
        cache_dir=tmp_path / "cache",
        use_gpu=False,
    )


@pytest.fixture
def sample_target() -> TargetProtein:
    return TargetProtein(
        identifier="P04629",
        name="TrkA (NGF receptor)",
        uniprot_id="P04629",
        sequence="MLRGGRRGQLGWHSWAAGPGSLLAWLI",  # truncated for testing
    )


@pytest.fixture
def sample_candidates() -> list[PeptideCandidate]:
    return [
        PeptideCandidate(sequence="ACDWKTHLYG", generation_method="pepmlm", generation_rank=1),
        PeptideCandidate(sequence="RLHFWKEPGM", generation_method="pepmlm", generation_rank=2),
        PeptideCandidate(sequence="QNWDCVRATK", generation_method="pepmlm", generation_rank=3),
    ]
