"""orthology.smk — OrthoFinder-based orthology clustering (Phase 2, supplementary).

OrthoFinder runs on ALL sequences across ALL families to produce de-novo orthogroups.
This is used to:
  (a) Validate that our enzyme_families.yaml accessions cluster as expected
  (b) Identify any additional family members missed in manual curation
  (c) Provide rooted gene trees for ancestral reconstruction

OrthoFinder output structure (created under results/orthologs/):
  OrthoFinder/Results_<date>/
    Orthogroups/Orthogroups.tsv
    Gene_Trees/
    Species_Tree/SpeciesTree_rooted.txt

NOTE: OrthoFinder discovers the output date-stamp directory name at runtime.
      The 'orthogroups_done' sentinel file is written by the shell block so
      downstream rules have a stable dependency target.
"""


rule run_orthofinder:
    """Run OrthoFinder on all downloaded sequences to cluster orthologs.

    Requires fetch_all_sequences checkpoint to have completed successfully.
    Input is the parent sequences directory; OrthoFinder discovers proteome
    files automatically via -f flag.
    """
    input:
        seqdir="results/sequences/",
    output:
        sentinel="results/orthologs/.orthofinder_done",
    params:
        inflation=config["orthology"]["inflation"],
        min_species=config["orthology"]["min_species"],
        outdir="results/orthologs/",
    conda:
        "../envs/bioinfo.yaml"
    log:
        "logs/orthology/orthofinder.log",
    threads: 16
    resources:
        mem_mb=32000,
        runtime=480,
    shell:
        """
        # OrthoFinder cannot write to an existing directory with prior results;
        # use a fresh subdirectory each run.
        mkdir -p {params.outdir}

        orthofinder \
            -f {input.seqdir} \
            -t {threads} \
            -a {threads} \
            -M msa \
            -I {params.inflation} \
            -o {params.outdir} \
            2> {log}

        # Write sentinel so downstream rules can depend on a stable path
        touch {output.sentinel}

        echo "OrthoFinder complete. Results in {params.outdir}" >> {log}
        """


rule extract_family_orthogroup:
    """Extract the orthogroup TSV rows for a single enzyme family.

    Uses the ANCHOR accessions from enzyme_families.yaml (e.g. AB114914 for
    nepenthesins) to identify the correct OrthoFinder orthogroup, then writes
    a filtered TSV.
    """
    input:
        sentinel="results/orthologs/.orthofinder_done",
        families_config="config/enzyme_families.yaml",
    output:
        tsv="results/orthologs/{family}_orthogroup.tsv",
    log:
        "logs/orthology/extract_{family}.log",
    shell:
        """
        python workflow/scripts/extract_orthogroup.py \
            --families-config {input.families_config} \
            --orthofinder-dir results/orthologs/ \
            --family {wildcards.family} \
            --output {output.tsv} \
            --verbose \
            2> {log}
        """
