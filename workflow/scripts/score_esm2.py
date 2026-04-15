#!/usr/bin/env python3
"""score_esm2 — ESM-2 masked marginal log-likelihood scoring of convergent substitutions.

For each convergent substitution identified in Phase 3, computes:
  LLR = log P(derived_aa | context) − log P(ancestral_aa | context)

where P(aa | context) is the ESM-2 model probability at the target position with all
other residues visible (masked marginal strategy, Meier et al. 2021, bioRxiv).

A positive LLR means the language model assigns higher probability to the derived
(carnivorous) amino acid than the ancestral residue — consistent with the derived
state being more tolerated in the sequence context. Combined with FoldX ΔΔG and
EVmutation ΔΔE, ESM-2 LLR forms the third independent axis of the convergent site
classification framework.

Reads:
  --convergent  results/convergence/{family}.convergent_sites.tsv
                Columns: alignment_position, ancestral_aa, derived_aa, species_list
  --sequences   results/sequences/{family}/all_sequences.fa
                All unique sequences in the family (used to derive sequence context)

Writes:
  --output  results/esm2/{family}.esm2_scores.tsv
            Columns: alignment_position, sequence_id, ancestral_aa, derived_aa,
                     log_p_ancestral, log_p_derived, llr, model

Citation:
  Lin Z, Akin H, Rao R, et al. (2023) Evolutionary-scale prediction of atomic-level
  protein structure with a language model. Science 379:1123–1130.
  Meier J, Rao R, Verkuil R, et al. (2021) Language models enable zero-shot prediction
  of the effects of mutations on protein function. NeurIPS 34.
"""

import logging
import sys
from pathlib import Path
from typing import Optional

import click
import numpy as np
import pandas as pd
import torch

logger = logging.getLogger(__name__)


def _load_esm_model(model_name: str, device: str):
    """Load ESM-2 model and batch converter. Downloads on first call."""
    try:
        import esm
    except ImportError:
        logger.error(
            "fair-esm not installed. Install with: pip install fair-esm"
        )
        sys.exit(1)

    logger.info("Loading ESM-2 model %s on device %s", model_name, device)
    model, alphabet = esm.pretrained.__dict__[model_name]()
    model = model.to(device)
    model.eval()
    batch_converter = alphabet.get_batch_converter()
    logger.info("ESM-2 model loaded (parameters: %dM)", sum(p.numel() for p in model.parameters()) // 1_000_000)
    return model, alphabet, batch_converter


def _compute_masked_marginals(
    model,
    alphabet,
    batch_converter,
    sequence: str,
    positions: list[int],
    device: str,
    batch_size: int = 32,
) -> dict[int, dict[str, float]]:
    """Return {position: {aa: log_prob}} for all positions in the list.

    Uses the masked marginal strategy: for each target position, mask it and
    compute log softmax over the vocabulary. This is the method from Meier et al.
    (2021, NeurIPS) and is the standard zero-shot fitness prediction approach.
    """
    mask_token = alphabet.mask_idx
    aa_tokens = {
        alphabet.get_idx(aa): aa
        for aa in "ACDEFGHIKLMNPQRSTVWY"
        if alphabet.get_idx(aa) != alphabet.unk_idx
    }

    results: dict[int, dict[str, float]] = {}

    # Batch positions to avoid OOM on large sequences
    for batch_start in range(0, len(positions), batch_size):
        batch_positions = positions[batch_start : batch_start + batch_size]
        seqs_masked = []
        for pos in batch_positions:
            seq_list = list(sequence)
            seq_list[pos] = "<mask>"
            seqs_masked.append(("seq", "".join(seq_list)))

        _, _, tokens = batch_converter(seqs_masked)
        tokens = tokens.to(device)

        with torch.no_grad():
            logits = model(tokens, repr_layers=[], return_contacts=False)["logits"]

        # logits shape: (batch, seq_len+2, vocab_size); +2 for BOS/EOS
        log_probs = torch.log_softmax(logits, dim=-1)

        for i, pos in enumerate(batch_positions):
            token_pos = pos + 1  # offset by BOS token
            pos_log_probs = log_probs[i, token_pos, :]
            results[pos] = {
                aa: pos_log_probs[tok_idx].item()
                for tok_idx, aa in aa_tokens.items()
            }

    return results


def _score_convergent_sites(
    convergent_df: pd.DataFrame,
    sequences: dict[str, str],
    model,
    alphabet,
    batch_converter,
    device: str,
    batch_size: int,
) -> pd.DataFrame:
    """Score each convergent substitution in every family sequence."""
    records = []

    for seq_id, sequence in sequences.items():
        logger.info("Scoring sequence %s (length %d)", seq_id, len(sequence))
        positions = convergent_df["sequence_position"].dropna().astype(int).tolist()
        positions = [p for p in positions if 0 <= p < len(sequence)]

        if not positions:
            logger.warning("No valid positions for sequence %s — skipping", seq_id)
            continue

        try:
            pos_scores = _compute_masked_marginals(
                model, alphabet, batch_converter, sequence, positions, device, batch_size
            )
        except RuntimeError as exc:
            logger.error("ESM-2 inference failed for %s: %s", seq_id, exc)
            continue

        for _, row in convergent_df.iterrows():
            pos = int(row["sequence_position"]) if pd.notna(row["sequence_position"]) else None
            if pos is None or pos not in pos_scores:
                continue
            anc = str(row["ancestral_aa"]).upper()
            der = str(row["derived_aa"]).upper()
            if anc not in pos_scores[pos] or der not in pos_scores[pos]:
                logger.warning(
                    "Amino acid %s or %s not in ESM vocabulary at position %d", anc, der, pos
                )
                continue
            log_p_anc = pos_scores[pos][anc]
            log_p_der = pos_scores[pos][der]
            records.append(
                {
                    "alignment_position": row["alignment_position"],
                    "sequence_position": pos,
                    "sequence_id": seq_id,
                    "ancestral_aa": anc,
                    "derived_aa": der,
                    "log_p_ancestral": log_p_anc,
                    "log_p_derived": log_p_der,
                    "llr": log_p_der - log_p_anc,
                }
            )

    return pd.DataFrame(records)


def _parse_fasta(fasta_path: Path) -> dict[str, str]:
    """Return {seq_id: sequence} from a FASTA file."""
    from Bio import SeqIO
    return {record.id: str(record.seq).replace("-", "").upper() for record in SeqIO.parse(fasta_path, "fasta")}


@click.command()
@click.option("--convergent", "-c", type=click.Path(exists=True), required=True,
              help="TSV of convergent sites from detect_convergence.py.")
@click.option("--sequences", "-s", type=click.Path(exists=True), required=True,
              help="FASTA with all family sequences (unaligned).")
@click.option("--output", "-o", type=click.Path(), required=True,
              help="Output TSV with ESM-2 LLR scores per site per sequence.")
@click.option("--model-name", default="esm2_t33_650M_UR50D",
              show_default=True,
              help="ESM-2 model name. Options: esm2_t6_8M_UR50D (fast/test), "
                   "esm2_t33_650M_UR50D (default), esm2_t36_3B_UR50D (slow/accurate).")
@click.option("--device", default="cuda", show_default=True,
              help="Torch device: 'cuda' or 'cpu'.")
@click.option("--batch-size", default=32, show_default=True,
              help="Positions per GPU batch for masked marginal inference.")
@click.option("--verbose", is_flag=True, help="Enable debug logging.")
def main(
    convergent: str,
    sequences: str,
    output: str,
    model_name: str,
    device: str,
    batch_size: int,
    verbose: bool,
) -> None:
    """Score convergent substitutions with ESM-2 masked marginal log-likelihood.

    Produces a log-likelihood ratio (LLR) for each convergent substitution in
    each family sequence. LLR > 0 means the language model favors the derived
    (carnivorous) amino acid over the ancestral state.
    """
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    convergent_path = Path(convergent)
    sequences_path = Path(sequences)
    output_path = Path(output)

    logger.info("Reading convergent sites from %s", convergent_path)
    convergent_df = pd.read_csv(convergent_path, sep="\t")
    required_cols = {"alignment_position", "ancestral_aa", "derived_aa", "sequence_position"}
    missing = required_cols - set(convergent_df.columns)
    if missing:
        logger.error("Missing columns in convergent TSV: %s", missing)
        sys.exit(1)
    logger.info("Loaded %d convergent sites", len(convergent_df))

    logger.info("Reading family sequences from %s", sequences_path)
    sequences_dict = _parse_fasta(sequences_path)
    logger.info("Loaded %d sequences", len(sequences_dict))

    if device == "cuda" and not torch.cuda.is_available():
        logger.warning("CUDA not available; falling back to CPU (slow)")
        device = "cpu"

    model, alphabet, batch_converter = _load_esm_model(model_name, device)

    logger.info("Starting masked marginal scoring")
    results_df = _score_convergent_sites(
        convergent_df, sequences_dict, model, alphabet, batch_converter, device, batch_size
    )

    if results_df.empty:
        logger.error("No scores produced — check input files and position mapping")
        sys.exit(1)

    results_df["model"] = model_name
    logger.info("Scored %d (site, sequence) pairs", len(results_df))

    # Summary stats for log
    pos_llr = results_df[results_df["llr"] > 0].shape[0]
    neg_llr = results_df[results_df["llr"] <= 0].shape[0]
    logger.info(
        "LLR summary: %d favorable (>0), %d unfavorable (≤0) out of %d",
        pos_llr, neg_llr, len(results_df),
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_path, sep="\t", index=False)
    logger.info("Wrote ESM-2 scores to %s", output_path)


if __name__ == "__main__":
    main()
