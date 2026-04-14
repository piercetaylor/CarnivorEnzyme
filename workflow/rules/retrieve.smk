"""retrieve.smk — Sequence download rules (Phase 1).

Two rules:
  1. fetch_all_sequences (checkpoint) — downloads all non-TODO accessions from
     enzyme_families.yaml via NCBI/UniProt into results/sequences/{family}/{species}.fa
  2. combine_family_sequences — concatenates per-species FASTAs into a single
     per-family combined.fa for MAFFT input.
"""

from pathlib import Path


def _get_species_fastas(wildcards):
    """Return list of per-species FASTA paths for a family after checkpoint resolves."""
    checkpoints.fetch_all_sequences.get()
    family_dir = Path(f"results/sequences/{wildcards.family}")
    # Exclude the combined.fa output itself if it already exists
    return sorted(
        p for p in family_dir.glob("*.fa")
        if p.name != "combined.fa"
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
        email=config["ncbi"]["email"],
        api_key=config["ncbi"]["api_key"],
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

    This is the direct input to MAFFT L-INS-i alignment.
    The lambda input triggers re-evaluation of the checkpoint DAG so that
    Snakemake discovers the actual FASTA files produced by fetch_all_sequences.
    """
    input:
        _get_species_fastas,
    output:
        combined="results/sequences/{family}/combined.fa",
    log:
        "logs/retrieve/combine_{family}.log",
    shell:
        """
        cat {input} > {output.combined} 2> {log}
        echo "Combined $(grep -c '^>' {output.combined}) sequences for family {wildcards.family}" \
            >> {log}
        """
