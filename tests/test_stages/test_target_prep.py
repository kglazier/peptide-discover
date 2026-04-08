"""Tests for target preparation stage."""

from pathlib import Path
from unittest.mock import patch

import pytest

from peptide_discover.models.target import TargetProtein
from peptide_discover.stages.target_prep import resolve_target


def test_resolve_unknown_identifier_raises():
    with pytest.raises(ValueError, match="Cannot resolve"):
        resolve_target("not_a_valid_id_at_all_xyz")


def test_resolve_from_file_fasta(tmp_path: Path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">test\nMKFLVLLFNILCSL\n")
    target = resolve_target(str(fasta))
    assert target.sequence == "MKFLVLLFNILCSL"
    assert target.structure_path is None


def test_resolve_from_file_unknown_format(tmp_path: Path):
    bad = tmp_path / "test.xyz"
    bad.write_text("junk")
    with pytest.raises(ValueError, match="Unsupported file format"):
        resolve_target(str(bad))


@pytest.mark.network
def test_resolve_uniprot_live():
    """Integration test — hits real UniProt + AlphaFold APIs."""
    target = resolve_target("P04629")  # TrkA, smaller protein
    assert len(target.sequence) > 100
    assert target.uniprot_id == "P04629"
    assert target.structure_path is not None
    assert target.structure_path.exists()
