"""Stage 2: Peptide candidate generation."""

import logging

from peptide_discover.models.peptide import PeptideCandidate
from peptide_discover.models.target import TargetProtein

logger = logging.getLogger(__name__)

PEPMLM_MODEL_ID = "TianlaiChen/PepMLM-650M"
PEPMLM_MAX_LENGTH = 550  # trained with max_length=552, minus BOS/EOS
STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")


def generate_pepmlm(
    target: TargetProtein,
    n: int = 100,
    peptide_length: int = 15,
    top_k: int = 3,
    device: str | None = None,
) -> list[PeptideCandidate]:
    """Generate peptide binder candidates using PepMLM.

    PepMLM is a fine-tuned ESM-2 (650M) model that generates linear peptide
    sequences conditioned on the target protein sequence. Single-pass masked
    prediction with top-k sampling for diversity.

    Args:
        target: Target protein with sequence.
        n: Number of candidates to generate.
        peptide_length: Length of peptides to generate (3-50).
        top_k: Top-k sampling parameter (higher = more diverse).
        device: Torch device string. Auto-detected if None.
    """
    try:
        import torch
        from torch.distributions import Categorical
        from transformers import AutoModelForMaskedLM, AutoTokenizer
    except ImportError:
        raise ImportError("PepMLM requires: pip install peptide-discover[pepmlm]")

    if not target.sequence:
        raise ValueError("Target protein must have a sequence for PepMLM generation.")

    if peptide_length < 3 or peptide_length > 50:
        raise ValueError("Peptide length must be between 3 and 50.")

    # Truncate target if needed to fit within model's max length
    max_protein_len = PEPMLM_MAX_LENGTH - peptide_length
    protein_seq = target.sequence[:max_protein_len]
    if len(protein_seq) < len(target.sequence):
        logger.warning(
            "Target sequence truncated from %d to %d residues to fit model max length.",
            len(target.sequence),
            len(protein_seq),
        )

    # Resolve device
    if device is None:
        from peptide_discover.utils.gpu import get_device
        device = get_device()

    logger.info("Loading PepMLM model on %s...", device)
    tokenizer = AutoTokenizer.from_pretrained(PEPMLM_MODEL_ID)
    model = AutoModelForMaskedLM.from_pretrained(PEPMLM_MODEL_ID).to(device)
    model.eval()

    # Build masked input: protein sequence + mask tokens for peptide
    masked_peptide = "<mask>" * peptide_length
    input_sequence = protein_seq + masked_peptide
    inputs = tokenizer(input_sequence, return_tensors="pt", truncation=True, max_length=552)
    inputs = {k: v.to(device) for k, v in inputs.items()}

    mask_token_id = tokenizer.mask_token_id
    mask_indices = (inputs["input_ids"] == mask_token_id).nonzero(as_tuple=True)[1]

    if len(mask_indices) == 0:
        raise RuntimeError("No mask tokens found in tokenized input.")

    # Generate candidates via repeated top-k sampling
    candidates = []
    seen_sequences = set()

    logger.info("Generating %d candidates (peptide_length=%d, top_k=%d)...", n, peptide_length, top_k)

    with torch.no_grad():
        logits = model(**inputs).logits
        logits_at_masks = logits[0, mask_indices]  # (peptide_length, vocab_size)
        top_k_logits, top_k_indices = logits_at_masks.topk(top_k, dim=-1)

    # Sample n candidates from the top-k distribution
    attempts = 0
    max_attempts = n * 10  # avoid infinite loop if diversity is exhausted

    with torch.no_grad():
        while len(candidates) < n and attempts < max_attempts:
            attempts += 1
            probabilities = torch.nn.functional.softmax(top_k_logits, dim=-1)
            sampled = Categorical(probabilities).sample()
            token_ids = top_k_indices.gather(-1, sampled.unsqueeze(-1)).squeeze(-1)
            sequence = tokenizer.decode(token_ids, skip_special_tokens=True).replace(" ", "")

            if sequence and sequence not in seen_sequences and all(
                aa in STANDARD_AAS for aa in sequence
            ):
                seen_sequences.add(sequence)
                candidates.append(
                    PeptideCandidate(
                        sequence=sequence,
                        generation_method="pepmlm",
                        generation_rank=len(candidates) + 1,
                    )
                )

    if len(candidates) < n:
        logger.warning(
            "Only generated %d unique candidates out of %d requested (top_k=%d may be too low).",
            len(candidates),
            n,
            top_k,
        )

    # Score candidates with pseudo-perplexity
    logger.info("Computing pseudo-perplexity scores for %d candidates...", len(candidates))
    candidates = _score_candidates(
        candidates, protein_seq, model, tokenizer, device
    )

    # Sort by perplexity (lower = better)
    candidates.sort(key=lambda c: c.generation_rank)

    logger.info("Generation complete. %d candidates.", len(candidates))
    return candidates


def _score_candidates(
    candidates: list[PeptideCandidate],
    protein_seq: str,
    model,
    tokenizer,
    device: str,
) -> list[PeptideCandidate]:
    """Score candidates using pseudo-perplexity.

    For each candidate, mask each peptide residue one at a time and measure
    how well the model reconstructs it. Lower perplexity = higher confidence.
    """
    import torch

    scored = []
    for candidate in candidates:
        full_seq = protein_seq + candidate.sequence
        inputs = tokenizer(full_seq, return_tensors="pt", truncation=True, max_length=552)
        inputs = {k: v.to(device) for k, v in inputs.items()}

        input_ids = inputs["input_ids"].clone()
        peptide_start = len(protein_seq) + 1  # +1 for BOS token
        peptide_end = peptide_start + len(candidate.sequence)

        log_likelihood = 0.0
        with torch.no_grad():
            for i in range(peptide_start, peptide_end):
                masked_ids = input_ids.clone()
                masked_ids[0, i] = tokenizer.mask_token_id
                logits = model(input_ids=masked_ids, attention_mask=inputs["attention_mask"]).logits
                probs = torch.nn.functional.softmax(logits[0, i], dim=-1)
                true_token = input_ids[0, i]
                log_likelihood += torch.log(probs[true_token]).item()

        # Perplexity = exp(-avg_log_likelihood)
        avg_ll = log_likelihood / len(candidate.sequence)
        perplexity = torch.exp(torch.tensor(-avg_ll)).item()

        scored.append(
            PeptideCandidate(
                sequence=candidate.sequence,
                generation_method=candidate.generation_method,
                generation_rank=0,  # will be reassigned after sorting
                length=candidate.length,
                perplexity=perplexity,
            )
        )

    # Sort by perplexity (lower = better) and assign ranks
    scored.sort(key=lambda c: c.perplexity or float("inf"))
    for i, c in enumerate(scored):
        c.generation_rank = i + 1

    return scored


def generate_rfdiffusion(target: TargetProtein, n: int = 100) -> list[PeptideCandidate]:
    """Generate peptide binder candidates using RFdiffusion + ProteinMPNN.

    Structure-based approach: generates backbone with RFdiffusion,
    then designs sequences with ProteinMPNN.
    """
    raise NotImplementedError("RFdiffusion generation not yet implemented.")
