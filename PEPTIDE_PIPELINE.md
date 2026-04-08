# AI Peptide Discovery Pipeline

## Goal

Build an integrated, consumer-hardware-friendly peptide discovery pipeline that orchestrates existing open-source tools into a single workflow. No new models — just the glue that's missing.

**Input:** A target protein (sequence or structure)
**Output:** Ranked novel peptide candidates with binding affinity scores AND safety/ADMET predictions, ready to hand to a CRO for synthesis and testing.

## Why This Doesn't Exist Yet

The individual tools are excellent. Nobody has wired them together:

- **ProteinDJ** (closest thing) requires HPC, focuses on protein binders, skips all property/safety prediction
- **PepMLM**, **Boltz-2**, **PeptiVerse**, **ToxinPred** etc. all work in isolation
- A researcher today has to manually run 5-6 tools, convert formats between them, and do their own ranking

## Pipeline Architecture

```
Target protein (sequence or PDB ID)
        |
        v
[1. Target Preparation]
    - Fetch structure from AlphaFold DB / PDB, or predict via ESMFold
    - Identify binding sites
        |
        v
[2. Peptide Generation]
    - PepMLM (sequence-only, linear peptides, consumer GPU)
    - RFdiffusion + ProteinMPNN (structure-based, more control)
    - Generate N candidates (default: 100-500)
        |
        v
[3. Binding Prediction]
    - Boltz-2 for structure prediction + affinity scoring
    - ~$0.01/prediction, single GPU, ~20 sec each
    - Filter: keep top K by predicted binding affinity
        |
        v
[4. Property / Safety Screening]
    - PeptiVerse: solubility, permeability, hemolysis, half-life, toxicity
    - ToxinPred 3.0: toxicity classification
    - B3Pred: blood-brain barrier penetration (critical for cognitive targets)
    - Composite safety score per candidate
        |
        v
[5. Ranking & Report]
    - Multi-objective ranking (binding affinity vs safety vs drug-likeness)
    - Comparison against known peptides for same target
    - Output: ranked list with all scores + CRO-ready specs
```

## Primary Research Tracks

### Track 1: Cognitive Enhancement
Target receptors:
- **TrkB** (BDNF pathway — neuroplasticity)
- **HGF/c-Met** (synaptogenesis — the Dihexa mechanism)
- **NGF receptors** (nerve growth factor)
- **AMPA receptor** modulators (learning/memory)

Success metric: candidates predicted to bind these targets with higher affinity than known peptides (Semax, Dihexa, etc.) while scoring lower on toxicity.

### Track 2: Muscle / Body Composition
Target receptors:
- **ActRIIB** (myostatin receptor — blocking = more muscle growth)
- **GHSR** (growth hormone secretagogue receptor — Ipamorelin target)
- **GHRH receptor** (CJC-1295 target)

Success metric: candidates predicted to bind with comparable or better affinity to known ligands, with improved stability/half-life.

### Track 3 (Future): Antimicrobial Peptides
Design novel AMPs targeting E. coli membrane proteins. **Direct integration with the AMR MIC prediction project** — design peptides here, predict their efficacy there.

This closed loop (design AMP → predict MIC against genomic resistance profiles) is potentially publishable and novel.

## Shared Infrastructure with AMR Project

Both projects use:
- **ESM-2** (but differently — see lesson learned below)
- **BioPython** for sequence handling
- The same architectural pattern: solved tool (AMRFinderPlus / PepMLM) + prediction layer on top

Shared code candidates:
- Sequence I/O and format conversion
- Results reporting

## Lesson Learned from AMR Project: ESM-2 Embedding Limitations

The AMR project tested ESM-2 mean-pooled embeddings as features for XGBoost MIC prediction. **Result: embeddings added noise that overwhelmed the signal.** Binary gene presence alone achieved 25/28 drugs above 90% essential agreement — better than the ESM-2-augmented model.

Root causes identified:
- Small training data (~500-650 samples per drug) — not enough for 32 compressed embedding features to help
- XGBoost can't learn from high-dimensional embeddings the way neural networks can
- Reference protein embeddings (not actual isolate sequences) miss the novel mutations ESM-2 is designed to capture

**How this applies to the peptide pipeline:**

The failure mode (mean-pooled embeddings → XGBoost → small N) does NOT apply here because the peptide tools use ESM-2 differently:
- **PepMLM**: fine-tuned ESM-2 doing masked sequence generation, not embedding extraction
- **Boltz-2 / PeptiVerse**: neural networks trained on millions of examples, not XGBoost on hundreds

But the broader principle carries forward:

> **Simple features that capture the right signal beat fancy embeddings that add noise.**

Our orchestration layer (step 5, ranking & report) should use **direct, interpretable scores** from each upstream tool (binding affinity number, toxicity flag, BBB yes/no) rather than attempting embedding-based meta-scoring or re-embedding candidates through additional models. Use the neural networks where they're already trained (inside PepMLM, Boltz-2, PeptiVerse). Keep the glue layer simple and transparent.

## Key Open-Source Dependencies

| Tool | Purpose | License | Hardware |
|------|---------|---------|----------|
| PepMLM (650M) | Peptide sequence generation | Open | Consumer GPU |
| RFdiffusion + ProteinMPNN | Structure-based peptide design | BSD / MIT | Consumer GPU (12GB+ VRAM) |
| Boltz-2 | Binding affinity + complex structure | MIT | Consumer GPU |
| ESMFold / AlphaFold DB | Target structure prediction | Open | Consumer GPU / API |
| PeptiVerse | Multi-property prediction (7 endpoints) | Open | Consumer GPU |
| ToxinPred 3.0 | Toxicity prediction | Open | CPU / Web |
| B3Pred | BBB penetration prediction | Open | CPU / Web |
| PepFuNN (Novo Nordisk) | Peptide analysis + characterization | Open | CPU |

## What "Done" Looks Like

### MVP
- CLI tool: `peptide-discover --target <uniprot_id_or_pdb> --track cognitive --candidates 100`
- Outputs a ranked CSV/JSON with: sequence, predicted binding affinity, toxicity score, BBB penetration, solubility, overall rank
- Runs on a single workstation with one GPU

### V2
- Web UI for non-technical users
- Batch mode for screening multiple targets
- Integration with AMR project for AMP design loop
- Support for non-natural amino acids (via PepINVENT / PepFoundry)
- Comparison dashboard against known peptides

## Realistic Expectations

- **High probability (~70-80%):** Find novel sequences that computationally bind targets
- **Moderate probability (~20-30%):** Candidates that validate in cell assays
- **Low probability (~5-10%):** Something meaningfully better than existing peptides
- **Lottery ticket (~0.1-1%):** A genuinely surprising discovery

The value is in the platform itself — it gets more powerful as models improve, and the AMP crossover with the AMR project is where the most novel research potential lives.
