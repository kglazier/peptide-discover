# Experimental Results

Real output from running the experiments in this folder. These are computational predictions, not wet-lab validated results. See the main [README.md](../README.md) for scope and limitations.

## Baseline: known peptides against real targets

From `experiments/baseline_known_peptides.py`. Docking known peptides against their real biological targets to calibrate what "good" scores look like.

| Peptide | Sequence | Target | Affinity (kcal/mol) |
|---|---|---|---|
| Ipamorelin | AIWFK | GHSR | **-8.33** |
| Magainin-2 | GIGKFLHS | LptD | -7.73 |
| Magainin-2 | GIGKFLHS | BamA | -7.02 |
| Selank | TKPRPGP | AGT | -6.07 |
| LL-37 (fragment) | LLGDFFRK | OmpA | -5.22 |
| Magainin-2 | GIGKFLHS | OmpA | -4.80 |
| Semax | MEHFPGP | POMC | -4.03 |
| Cecropin-A (fragment) | KWKLFKKI | OmpA | 0.00 (dock failed) |

**Interpretation:** Ipamorelin binding GHSR at -8.33 is a real drug hitting its real receptor. This calibrates the pipeline — scores in the -8 range are meaningful.

## AMP design for E. coli membrane targets

From `experiments/amp_mdr_ecoli.py`. 30 candidates generated per target at 8 residues, docked with AutoDock Vina at exhaustiveness=4.

### Summary by target

| Target | UniProt | Size (aa) | Best binder | Score | Success rate | Magainin-2 ref |
|---|---|---|---|---|---|---|
| **LptD** | P31554 | 784 | WSAGLSLT | **-8.65** | 30/30 | -7.73 |
| **BamA** | P0A940 | 810 | VVPSAMSA | -7.82 | 30/30 | -7.02 |
| **OmpA** | P0A910 | 346 | DNDAAKKN | -5.83 | 16/30 | -4.80 |

**Key finding:** Generated candidates score comparable-or-better than natural Magainin-2 against BamA and LptD. OmpA is a difficult target for both natural and generated peptides — likely lacks a clean binding pocket.

### Top 5 candidates per target

#### LptD (best target)

| Rank | Sequence | Affinity | Toxicity | Solubility | Notes |
|---|---|---|---|---|---|
| 1 | WSAGLSLT | -8.65 | 0.00 | 0.41 | Strongest binder overall |
| 2 | DANSAKLA | -8.28 | 0.00 | 0.54 | Clean profile |
| 3 | WSAGSSLA | -8.10 | 0.00 | 0.44 | |
| 4 | SRYGGSKS | -8.07 | 0.10 | **0.90** | Best solubility in top 5 |
| 5 | GSASASLS | -8.01 | 0.00 | 0.44 | |

#### BamA

| Rank | Sequence | Affinity | Toxicity | Solubility | Notes |
|---|---|---|---|---|---|
| 1 | VVPSAMSA | -7.82 | 0.00 | 0.13 | Strong binding, poor solubility |
| 2 | TAAEEGGA | -7.58 | 0.00 | 0.75 | **Best overall profile for BamA** |
| 3 | AVESLGLS | -7.56 | 0.00 | 0.47 | Lowest perplexity (model confident) |
| 4 | SDEINGLA | -7.50 | 0.00 | 0.73 | Clean profile |
| 5 | SVEEALLG | -7.49 | 0.00 | 0.62 | |

#### OmpA (hard target)

| Rank | Sequence | Affinity | Toxicity | Solubility | Notes |
|---|---|---|---|---|---|
| 1 | DNDAAKKN | -5.83 | 0.00 | 0.78 | Best of a weak field |
| 2 | DVVAQADA | -5.77 | 0.00 | 0.65 | |
| 3 | DEDAQKGK | -5.73 | 0.00 | 0.92 | High solubility |
| 4 | GEPADALA | -5.58 | 0.00 | 0.70 | |
| 5 | DANAQQGK | -5.51 | 0.00 | 0.73 | |

**Note:** OmpA docking had 14/30 failures (affinity=0.0). The structure appears poorly suited for small-molecule-like docking.

## Resistance context (from AMRCast)

Analysis of 10,654 E. coli genomes from the NARMS antibiogram dataset (via [AMRCast](https://github.com/kglazier/amrcast)):

- **53.5%** show resistance to at least one antibiotic
- **21.3%** are resistant to 5+ antibiotics (MDR)
- **Top resistance:** ampicillin (35.5%), tetracycline (29.5%), ceftriaxone (14.4%)

This motivates the search for AMPs that target essential membrane proteins — resistance mechanisms that work against traditional antibiotics (beta-lactamases, efflux pumps) don't protect against membrane disruption.

## How to reproduce

```bash
# Baseline (~15 min)
python experiments/baseline_known_peptides.py

# Full AMP design (~90 min per target)
python experiments/amp_mdr_ecoli.py -n 30 -l 8 -e 4 --targets BamA,LptD
```

Results are written to `results/` (gitignored).

## Caveats

- **Docking scores aren't truth.** Vina approximates free energy but isn't perfect. Two candidates with similar scores may behave very differently in practice.
- **No wet lab validation.** These are computational predictions only. Actual antimicrobial activity requires MIC assays.
- **8-mer peptides.** Shorter peptides dock faster but may be less specific. Longer peptides (12+) showed much higher Vina failure rates on consumer hardware.
- **Safety screening is property-based, not ML-based.** Catches obvious red flags but won't predict subtle toxicity.
