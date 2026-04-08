"""Predefined research track configurations."""

from dataclasses import dataclass, field


@dataclass
class ResearchTrack:
    name: str
    description: str
    targets: dict[str, str]  # friendly name -> UniProt ID
    require_bbb: bool = False
    reference_peptides: list[str] = field(default_factory=list)


TRACKS: dict[str, ResearchTrack] = {
    "cognitive": ResearchTrack(
        name="cognitive",
        description="Cognitive enhancement — neuroplasticity and synaptogenesis targets",
        targets={
            "TrkB": "Q16620",       # BDNF receptor
            "c-Met": "P08581",      # HGF receptor (Dihexa target)
            "TrkA": "P04629",       # NGF receptor
            "GluA1": "P42261",      # AMPA receptor subunit
        },
        require_bbb=True,
        reference_peptides=["Semax", "Selank", "Dihexa"],
    ),
    "muscle": ResearchTrack(
        name="muscle",
        description="Muscle growth and body composition targets",
        targets={
            "ActRIIB": "Q13705",    # Myostatin / activin receptor
            "GHSR": "Q92847",       # Growth hormone secretagogue receptor (Ipamorelin)
            "GHRHR": "Q02643",      # GHRH receptor (CJC-1295)
        },
        require_bbb=False,
        reference_peptides=["Ipamorelin", "CJC-1295", "Follistatin-344"],
    ),
    "amp": ResearchTrack(
        name="amp",
        description="Antimicrobial peptides — E. coli membrane targets",
        targets={
            "OmpA": "P0A910",       # Outer membrane protein A
            "BamA": "P0A940",       # Outer membrane assembly factor
            "LptD": "P31554",       # LPS transport protein
        },
        require_bbb=False,
        reference_peptides=["Magainin-2", "Cecropin-A", "LL-37"],
    ),
}


def get_track(name: str) -> ResearchTrack:
    """Get a research track by name."""
    if name not in TRACKS:
        available = ", ".join(TRACKS.keys())
        raise ValueError(f"Unknown track '{name}'. Available: {available}")
    return TRACKS[name]
