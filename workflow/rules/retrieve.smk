"""retrieve.smk — Sequence download rules (Phase 1).

Two rules:
  1. fetch_all_sequences (checkpoint) — downloads all non-TODO accessions from
     enzyme_families.yaml via NCBI/UniProt into results/sequences/{family}/{species}.fa
  2. combine_family_sequences — concatenates per-species FASTAs into a single
     per-family FASTA for MAFFT input, written to results/family_fasta/.

Path note (fixed 2026-09-02, audit/16): combine_family_sequences MUST NOT write inside
results/sequences/, because that directory is the declared output of the
fetch_all_sequences checkpoint. Snakemake rejects a DAG in which one rule's output is a
child of another's with ChildIOException, which blocked every dry run past Phase 1.
"""

import os
from pathlib import Path


def _get_species_fastas(wildcards):
    """Return list of per-species FASTA paths for a family after checkpoint resolves."""
    checkpoints.fetch_all_sequences.get()
    family_dir = Path(f"results/sequences/{wildcards.family}")
    # Defensive: combine_family_sequences now writes to results/family_fasta/, so a
    # combined.fa here can only be a leftover from a pre-2026-09-02 run. Excluding it
    # keeps a stale artefact from being concatenated into its own successor.
    return sorted(
        p for p in family_dir.glob("*.fa")
        if p.name not in {"combined.fa", "all_sequences.fa"}
    )


checkpoint fetch_all_sequences:
    """Download all protein sequences listed in enzyme_families.yaml.

    Output is a directory; exact files are not known until runtime because
    they depend on which YAML accessions resolve successfully (TODO entries
    are skipped; failed fetches are logged but do not abort if only partial
    failures occur in non-required species).
    """
    input:
        families_config="config/enzyme_families.yaml",
    output:
        seqdir=directory("results/sequences/"),
    params:
        email=config["ncbi"]["email"] or os.environ.get("NCBI_EMAIL", ""),
        api_key=config["ncbi"]["api_key"] or os.environ.get("NCBI_API_KEY", ""),
    conda:
        "../envs/bioinfo.yaml"
    log:
        "logs/retrieve/fetch_sequences.log",
    threads: 1
    shell:
        """
        python workflow/scripts/fetch_sequences.py \
            --families-config {input.families_config} \
            --output-dir {output.seqdir} \
            --email "{params.email}" \
            --api-key "{params.api_key}" \
            --verbose \
            2> {log}
        """


rule combine_family_sequences:
    """Concatenate all per-species FASTAs for one family into a single file.

    This is the direct input to MAFFT L-INS-i alignment (alignment.smk) and to ESM-2
    scoring (esm2.smk). Output lives OUTSIDE results/sequences/ — see the path note in
    this file's module docstring.
    The lambda input triggers re-evaluation of the checkpoint DAG so that
    Snakemake discovers the actual FASTA files produced by fetch_all_sequences.
    """
    input:
        _get_species_fastas,
    output:
        combined="results/family_fasta/{family}.combined.fa",
    log:
        "logs/retrieve/combine_{family}.log",
    shell:
        """
        cat {input} > {output.combined} 2> {log}
        echo "Combined $(grep -c '^>' {output.combined}) sequences for family {wildcards.family}" \
            >> {log}
        """
