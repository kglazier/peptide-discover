"""Central configuration for peptide-discover."""

from pathlib import Path

from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    model_config = {"env_prefix": "PEPDISC_"}

    # Paths
    data_dir: Path = Path("data")
    results_dir: Path = Path("results")
    cache_dir: Path = Path(".cache/peptide_discover")

    # GPU
    use_gpu: bool = True
    gpu_device: int = 0
    max_vram_gb: float = 12.0

    # Generation
    default_candidates: int = 100
    default_method: str = "pepmlm"

    # Binding (Boltz-2)
    binding_top_k: int = 50

    # Screening toggles
    enable_toxinpred: bool = True
    enable_b3pred: bool = True
    enable_peptiverse: bool = True


def get_settings() -> Settings:
    """Load settings from environment."""
    return Settings()
