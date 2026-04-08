# peptide-discover

AI peptide discovery pipeline — orchestrating generation, binding prediction, and safety screening.

## Install

```bash
pip install -e ".[dev]"
```

## Usage

```bash
peptide-discover --version
peptide-discover run pipeline --target Q16620 --track cognitive --candidates 100
```

See [PEPTIDE_PIPELINE.md](PEPTIDE_PIPELINE.md) for full architecture and design decisions.
