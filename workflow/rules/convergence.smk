rule detect_convergence:
    """Detect convergent amino acid substitutions from ancestral reconstruction.

    Consumes the tree AND the .state file produced by the SAME IQ-TREE --ancestral
    invocation (rule ancestral_reconstruction in phylogeny.smk). Passing
    {family}.rooted.treefile here instead would reintroduce the node-correspondence bug
    documented in audit/11_ancestral_reconstruction_architecture_fix.md: IQ-TREE's NodeN
    labels describe the tree IT wrote, and rerooting moves them onto other clades.
    """
    input:
        alignment="results/alignments/{family}.trimmed.fa",
        tree="results/phylogenies/{family}.asr.treefile",
        state="results/phylogenies/{family}.asr.state",
        species_yaml="config/species.yaml",
    output:
        tsv="results/convergence/{family}.convergent_sites.tsv",
    params:
        threshold=config["convergence"]["posterior_threshold"],
        min_lineages=config["convergence"]["min_lineages"],
    conda:
        "../envs/bioinfo.yaml"
    log:
        "logs/convergence/{family}.log",
    threads: 1
    shell:
        """
        python workflow/scripts/detect_convergence.py \
            --alignment {input.alignment} \
            --tree {input.tree} \
            --state {input.state} \
            --species-yaml {input.species_yaml} \
            --output {output.tsv} \
            --family {wildcards.family} \
            --threshold {params.threshold} \
            --min-lineages {params.min_lineages} \
            --verbose \
            2> {log}
        """
