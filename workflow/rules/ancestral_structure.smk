"""ancestral_structure.smk — Phase 4A: Ancestral sequence reconstruction + structure prediction.

Reconstructs the last common ancestor (MRCA) of all carnivorous plant taxa for each
enzyme family using IQ-TREE marginal ancestral state reconstruction (.state file from
Phase 2). Predicts the ancestral structure with AlphaFold3 / Chai-1, then compares
the ancestral structure to each modern ortholog to quantify structural divergence at
convergent sites.

Scientific rationale:
  Comparing the ancestral (pre-carnivory) enzyme structure to modern carnivorous
  plant structures identifies which structural features were acquired during the
  independent carnivorous transitions. If convergent substitutions produce coordinated
  active-site rearrangements, the RMSD restricted to convergent positions (RMSD_conv)
  will significantly exceed background RMSD — implicating these sites in functional
  restructuring rather than neutral drift.

Rules:
  1. extract_ancestor    — Parse IQ-TREE .state → ancestral FASTA
  2. predict_ancestor    — AF3/Chai-1 structure prediction for ancestral sequence
  3. compare_structures  — TM-align + convergent-site RMSD vs all modern structures

Citations:
  Pupko T et al. (2000) Mol Biol Evol 17:890.
  Hanson-Smith V, Johnson A. (2016) Low-complexity ancestral sequence reconstruction.
  Elife 5:e12094. [PASR framework; our approach is similar]
"""


rule extract_ancestor:
    """Parse IQ-TREE .state file and write MRCA ancestral sequence FASTA.

    Identifies the MRCA of carnivorous taxa in the rooted gene tree, builds the
    MAP (maximum a posteriori) sequence, and writes it for structure prediction.
    Sites with posterior probability < config.ancestral.mrca_posterior_threshold
    are masked as 'X' and excluded from downstream FoldX analysis.
    """
    input:
        state="results/phylogenies/{family}.state",
        tree="results/phylogenies/{family}.rooted.treefile",
    output:
        fasta="results/ancestral/{family}.mrca_ancestor.fa",
        stats="results/ancestral/{family}.mrca_stats.tsv",
    params:
        threshold=config["ancestral"]["mrca_posterior_threshold"],
    conda:
        "../envs/bioinfo.yaml"
    log:
        "logs/ancestral/{family}.extract.log",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=15,
    shell:
        """
        python workflow/scripts/extract_ancestor.py \
            --state-file          {input.state} \
            --tree-file           {input.tree} \
            --family              {wildcards.family} \
            --output              {output.fasta} \
            --output-stats        {output.stats} \
            --posterior-threshold {params.threshold} \
            --verbose \
            2> {log}
        """


rule predict_ancestor_chai1:
    """Predict ancestral enzyme structure with Chai-1 (fallback predictor).

    Runs when AF3 is unavailable or during development iteration. Outputs
    mmCIF format; assess_structure.py extracts pLDDT and converts to PDB.
    """
    input:
        fasta="results/ancestral/{family}.mrca_ancestor.fa",
    output:
        cif="results/ancestral/{family}.mrca_ancestor.cif",
        pdb="results/ancestral/{family}.mrca_ancestor.pdb",
        plddt="results/ancestral/{family}.mrca_ancestor.plddt.tsv",
    conda:
        "../envs/structure.yaml"
    log:
        "logs/ancestral/{family}.predict.log",
    threads: 4
    resources:
        mem_mb=32000,
        runtime=120,
        gpu=1,
    shell:
        """
        python workflow/scripts/predict_chai1.py \
            --input   {input.fasta} \
            --output  {output.cif} \
            --n-seeds 1 \
            --verbose \
            2> {log}

        python workflow/scripts/assess_structure.py \
            --input    {output.cif} \
            --output   {output.pdb} \
            --plddt    {output.plddt} \
            --verbose \
            2>> {log}
        """


rule compare_ancestral_modern:
    """Compare ancestral and modern carnivorous plant structures per family.

    Runs TM-align for global structural comparison and computes Cα RMSD
    restricted to convergent substitution positions. The question: are convergent
    sites structurally displaced between ancestor and modern, consistent with
    functional restructuring?
    """
    input:
        ancestral_pdb="results/ancestral/{family}.mrca_ancestor.pdb",
        modern_dir="results/structures/{family}/",
        convergent="results/convergence/{family}.convergent_sites.tsv",
    output:
        tsv="results/ancestral/{family}.structural_comparison.tsv",
    params:
        active_site_radius=config["ancestral"]["active_site_radius_angstrom"],
    conda:
        "../envs/structure.yaml"
    log:
        "logs/ancestral/{family}.compare.log",
    threads: 2
    resources:
        mem_mb=8000,
        runtime=60,
    shell:
        """
        python workflow/scripts/compare_ancestor_modern.py \
            --ancestral-pdb      {input.ancestral_pdb} \
            --modern-dir         {input.modern_dir} \
            --convergent         {input.convergent} \
            --output             {output.tsv} \
            --active-site-radius {params.active_site_radius} \
            --verbose \
            2> {log}
        """
