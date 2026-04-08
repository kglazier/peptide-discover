"""GPU detection and management."""


def detect_gpu() -> dict:
    """Detect available GPU and return info."""
    try:
        import torch

        if torch.cuda.is_available():
            device = torch.cuda.get_device_name(0)
            vram = torch.cuda.get_device_properties(0).total_memory / (1024**3)
            return {
                "available": True,
                "device": device,
                "vram_gb": round(vram, 1),
                "cuda_version": torch.version.cuda or "unknown",
            }
    except ImportError:
        pass

    return {"available": False, "device": None, "vram_gb": 0, "cuda_version": None}


def get_device() -> str:
    """Return the best available torch device string."""
    info = detect_gpu()
    return "cuda:0" if info["available"] else "cpu"


def check_vram(required_gb: float) -> bool:
    """Check if enough VRAM is available."""
    info = detect_gpu()
    return info["vram_gb"] >= required_gb
