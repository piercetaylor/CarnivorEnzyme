#!/usr/bin/env python3
"""fetch_sequences — Download protein sequences from NCBI and UniProt.

Reads config/enzyme_families.yaml, fetches each verified accession via the
Entrez API (NCBI) or UniProt REST, validates sequence lengths against expected
family ranges, and writes per-family per-species FASTA files to:
    results/sequences/{family}/{species}.fa

Skip logic:
  - Entries starting with "TODO" are silently skipped.
  - Sequences shorter than min_length_fraction * median for that family
    are excluded with a WARNING log line.

Output files written:
  - {output_dir}/{family}/{species}.fa  (one FASTA per species, ≥1 record)
"""

import logging
import re
import sys
import time
from io import StringIO
from pathlib import Path
from statistics import median
from typing import Optional

import click
import requests
import yaml
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

# ── Accession format detection ────────────────────────────────────────────────
# UniProt Swiss-Prot: [A-Z][0-9][A-Z0-9]{3}[0-9]  (6 chars)
# UniProt TrEMBL:     [A-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9] (10 chars)
_UNIPROT_RE = re.compile(
    r"^[A-Z][0-9][A-Z0-9]{3}[0-9]$"
    r"|^[A-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]$"
)

# UniProt FASTA endpoint (REST API v2)
_UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"


def _is_uniprot(accession: str) -> bool:
    """Return True if accession matches UniProt Swiss-Prot or TrEMBL format."""
    return bool(_UNIPROT_RE.match(accession.split(".")[0]))


# ── Fetch helpers ─────────────────────────────────────────────────────────────

def _fetch_ncbi(accession: str, email: str, api_key: Optional[str]) -> Optional[SeqRecord]:
    """Fetch one protein sequence from NCBI Entrez."""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    try:
        handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return record
    except Exception as exc:
        logger.warning("NCBI fetch failed for %s: %s", accession, exc)
        return None


def _fetch_uniprot(accession: str) -> Optional[SeqRecord]:
    """Fetch one protein sequence from UniProt REST API."""
    base_acc = accession.split(".")[0]  # strip version suffix
    url = _UNIPROT_URL.format(acc=base_acc)
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        record = SeqIO.read(StringIO(resp.text), "fasta")
        return record
    except Exception as exc:
        logger.warning("UniProt fetch failed for %s: %s", accession, exc)
        return None


def _fetch(
    accession: str,
    email: str,
    api_key: Optional[str],
    rate_delay: float,
) -> Optional[SeqRecord]:
    """Dispatch fetch to NCBI or UniProt; honour rate limit delay."""
    time.sleep(rate_delay)
    if _is_uniprot(accession):
        logger.debug("→ UniProt: %s", accession)
        return _fetch_uniprot(accession)
    else:
        logger.debug("→ NCBI: %s", accession)
        return _fetch_ncbi(accession, email, api_key)


# ── Validation helpers ────────────────────────────────────────────────────────

def _is_todo(value: str) -> bool:
    """Return True if accession is a placeholder (starts with TODO or is empty)."""
    return not value or str(value).strip().upper().startswith("TODO")


# ── Main ──────────────────────────────────────────────────────────────────────

@click.command()
@click.option(
    "--families-config", "-f",
    type=click.Path(exists=True),
    required=True,
    help="Path to enzyme_families.yaml.",
)
@click.option(
    "--output-dir", "-o",
    type=click.Path(),
    required=True,
    help="Root directory for output FASTAs (e.g. results/sequences/).",
)
@click.option(
    "--email", "-e",
    type=str,
    required=True,
    help="NCBI Entrez email address.",
)
@click.option(
    "--api-key", "-k",
    type=str,
    default="",
    help="NCBI API key (enables 10 req/s; 3 req/s without key).",
)
@click.option(
    "--families",
    type=str,
    default="",
    help="Comma-separated list of tier1 family keys to fetch (default: all).",
)
@click.option(
    "--min-length-fraction",
    type=float,
    default=0.8,
    show_default=True,
    help="Minimum sequence length as fraction of family median length.",
)
@click.option("--verbose", is_flag=True, help="Enable DEBUG logging.")
def main(
    families_config: str,
    output_dir: str,
    email: str,
    api_key: str,
    families: str,
    min_length_fraction: float,
    verbose: bool,
) -> None:
    """Download protein sequences for all enzyme families in enzyme_families.yaml."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    config_path = Path(families_config)
    output_root = Path(output_dir)

    with config_path.open(encoding="utf-8") as fh:
        config = yaml.safe_load(fh)

    tier1: dict = config.get("tier1", {})
    family_filter: set[str] = {f.strip() for f in families.split(",") if f.strip()}
    if family_filter:
        logger.info("Restricting to families: %s", sorted(family_filter))

    # Rate: 3 req/s without key, 10/s with key
    rate_delay = 0.10 if api_key else 0.34

    n_fetched = 0
    n_skipped = 0
    n_failed = 0

    for family_key, family_data in tier1.items():
        if family_filter and family_key not in family_filter:
            continue

        exp_min, exp_max = family_data.get("expected_length_aa", [0, 99999])
        species_map: dict = family_data.get("accessions", {})
        family_dir = output_root / family_key
        family_dir.mkdir(parents=True, exist_ok=True)

        logger.info("━━━ Family: %s (expected %d–%d aa) ━━━", family_key, exp_min, exp_max)

        # ── Pass 1: fetch all non-TODO accessions ─────────────────────────────
        fetched: dict[str, list[SeqRecord]] = {}

        for species, acc_list in species_map.items():
            if not isinstance(acc_list, list):
                continue
            species_records: list[SeqRecord] = []

            for acc in acc_list:
                acc = str(acc).strip()
                if _is_todo(acc):
                    logger.debug("SKIP TODO  %s / %s", family_key, species)
                    n_skipped += 1
                    continue

                record = _fetch(acc, email, api_key or None, rate_delay)
                if record is None:
                    logger.error("FAILED     %s / %s / %s", family_key, species, acc)
                    n_failed += 1
                    continue

                seq_len = len(record.seq)
                if seq_len < exp_min * 0.5:
                    logger.warning(
                        "VERY SHORT %s: %d aa (expected %d–%d) — included but flagged",
                        acc, seq_len, exp_min, exp_max,
                    )
                elif seq_len < exp_min:
                    logger.warning(
                        "SHORT      %s: %d aa (expected >=%d aa) — flagged as partial",
                        acc, seq_len, exp_min,
                    )

                # Normalise record ID and description for downstream parsing
                record.id = f"{species}|{acc}"
                record.name = acc
                record.description = (
                    f"family={family_key} species={species} acc={acc} len={seq_len}"
                )
                species_records.append(record)
                n_fetched += 1
                logger.info("  FETCHED  %s | %s (%d aa)", species, acc, seq_len)

            if species_records:
                fetched[species] = species_records

        # ── Pass 2: median-length filter ─────────────────────────────────────
        all_lengths = [len(r.seq) for recs in fetched.values() for r in recs]
        if not all_lengths:
            logger.warning("No sequences fetched for family %s — check YAML TODOs", family_key)
            continue

        median_len = median(all_lengths)
        cutoff = min_length_fraction * median_len
        logger.info(
            "Family %s: %d sequences, median %d aa, length cutoff %.0f aa",
            family_key, len(all_lengths), round(median_len), cutoff,
        )

        for species, records in fetched.items():
            kept: list[SeqRecord] = []
            for rec in records:
                seq_len = len(rec.seq)
                if seq_len >= cutoff:
                    kept.append(rec)
                else:
                    logger.warning(
                        "FILTERED   %s (%d aa) < %.0f%% of median (%.0f aa)",
                        rec.id, seq_len, min_length_fraction * 100, median_len,
                    )
                    n_skipped += 1

            if not kept:
                continue
            out_path = family_dir / f"{species}.fa"
            SeqIO.write(kept, str(out_path), "fasta")
            logger.info("  WROTE    %s → %s (%d record(s))", species, out_path, len(kept))

    logger.info(
        "Done. fetched=%d  skipped/filtered=%d  failed=%d",
        n_fetched, n_skipped, n_failed,
    )
    if n_failed > 0:
        logger.error(
            "%d accession(s) failed — re-run with --verbose to debug; "
            "check NCBI/UniProt availability for those IDs",
            n_failed,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
