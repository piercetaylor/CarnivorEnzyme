"""phylogeny.smk — Maximum-likelihood phylogenetic inference (Phase 2).

One rule per family: infer_tree
  - IQ-TREE 3 with ModelFinder Plus (-m MFP) for automatic model selection
  - UFBoot2 (-bb 1000) for ultrafast bootstrap approximation
  - SH-aLRT (-alrt 1000) for additional branch support metric
  - -ancestral for marginal ancestral state reconstruction (required by Phase 3)
  - -seed 42 for reproducibility

Input:  results/alignments/{family}.trimmed.fa
Output:
  results/phylogenies/{family}.treefile  — ML tree with branch supports
  results/phylogenies/{family}.state     — marginal ancestral states (used by detect_convergence.py)

NOTE on binary name:
  IQ-TREE 3 installs as 'iqtree3' in the conda package. If the cluster has
  an older module loaded as 'iqtree', substitute accordingly. The bioinfo.yaml
  conda environment pins iqtree>=3.0.0 which provides the iqtree3 binary.

Tool justification (CLAUDE.md §1):
  IQ-TREE 3 (Minh et al. 2025, EcoEvoRxiv): improved site concordance factors
  that reduce homoplasy effects in ancestral reconstruction — directly relevant
  to Phase 3 convergence detection. piqtree Python bindings enable integration
  with detect_convergence.py.
"""


rule infer_tree:
    """Infer ML phylogeny and ancestral states with IQ-TREE 3.

    Key flags:
      -m MFP        ModelFinder Plus: automatically selects best-fit model
      -bb 1000      UFBoot2 ultrafast bootstrap (1000 replicates)
      -alrt 1000    SH-aLRT branch test (report alongside UFBoot2)
      -seed 42      Reproducible random seed
      -ancestral    Marginal ancestral state reconstruction (writes .state file)
      -redo         Overwrite previous run in same prefix directory
      -pre          Output prefix (all IQ-TREE files share this stem)

    UFBoot2 interpretation: values ≥95 indicate strong support (NOT comparable
    to standard bootstrap ≥70 cutoff). Use alrt ≥80 AND ufboot ≥95 together.
    """
    input:
        alignment="results/alignments/{family}.trimmed.fa",
    output:
        treefile="results/phylogenies/{family}.treefile",
        state="results/phylogenies/{family}.state",
        iqlog="results/phylogenies/{family}.iqtree",
    params:
        bootstrap=config["phylogeny"]["bootstrap"],
        model=config["phylogeny"]["model"],
        prefix="results/phylogenies/{family}",
    conda:
        "../envs/bioinfo.yaml"
    log:
        "logs/phylogeny/{family}.log",
    threads: 8
    resources:
        mem_mb=16000,
        runtime=240,
    shell:
        """
        mkdir -p results/phylogenies/

        iqtree3 \
            -s {input.alignment} \
            -m {params.model} \
            -bb {params.bootstrap} \
            -alrt 1000 \
            -nt {threads} \
            -seed 42 \
            -ancestral \
            -redo \
            -pre {params.prefix} \
            2> {log}

        # Confirm required output files exist
        test -f {output.treefile} || \
            {{ echo "ERROR: treefile not produced by IQ-TREE" >> {log}; exit 1; }}
        test -f {output.state} || \
            {{ echo "ERROR: .state file not produced — check -ancestral flag" >> {log}; exit 1; }}

        echo "IQ-TREE complete for {wildcards.family}" >> {log}
        grep "Best-fit model:" {output.iqlog} >> {log} 2>/dev/null || true
        """


rule root_tree:
    """Root the IQ-TREE unrooted tree using outgroup sequences.

    Uses the species listed as outgroup_species in species.yaml to identify
    appropriate outgroup taxa. Writes a rooted Newick tree used in Phase 3
    ancestral reconstruction.

    Requires piqtree (IQ-TREE 3 Python bindings) or ete3 for tree manipulation.
    """
    input:
        treefile="results/phylogenies/{family}.treefile",
        species_config="config/species.yaml",
        families_config="config/enzyme_families.yaml",
    output:
        rooted="results/phylogenies/{family}.rooted.treefile",
    conda:
        "../envs/bioinfo.yaml"
    log:
        "logs/phylogeny/{family}.root.log",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10,
    shell:
        """
        python workflow/scripts/root_tree.py \
            --treefile {input.treefile} \
            --species-config {input.species_config} \
            --families-config {input.families_config} \
            --family {wildcards.family} \
            --output {output.rooted} \
            --verbose \
            2> {log}
        """
