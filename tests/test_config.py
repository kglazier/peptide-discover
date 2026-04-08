"""Test configuration loading."""

from peptide_discover.config.settings import Settings, get_settings
from peptide_discover.config.tracks import TRACKS, get_track


def test_default_settings():
    s = Settings(use_gpu=False)
    assert s.default_candidates == 100
    assert s.default_method == "pepmlm"
    assert s.binding_top_k == 50


def test_get_settings():
    s = get_settings()
    assert isinstance(s, Settings)


def test_tracks_exist():
    assert "cognitive" in TRACKS
    assert "muscle" in TRACKS
    assert "amp" in TRACKS


def test_cognitive_track():
    t = get_track("cognitive")
    assert t.require_bbb is True
    assert "TrkB" in t.targets


def test_muscle_track():
    t = get_track("muscle")
    assert t.require_bbb is False
    assert "GHSR" in t.targets


def test_unknown_track_raises():
    import pytest

    with pytest.raises(ValueError, match="Unknown track"):
        get_track("nonexistent")
