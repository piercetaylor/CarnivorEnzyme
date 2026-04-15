"""fep.smk — Phase 5D: Alchemical FEP for substrate binding ΔΔG.

Replaces the fep stub. GROMACS + pmx compute ΔΔG_bind for the top docking
substrates at each convergent substitution site, using the thermodynamic cycle:

  ΔΔG_bind = ΔG_holo(ancestral→derived) − ΔG_apo(ancestral→derived)

A negative ΔΔG_bind means the derived (carnivorous) amino acid tightens substrate
binding — mechanistic evidence for adaptive substrate specificity change.

Protocol (13 λ-windows × 5 ns × 3 replicates; MBAR estimator):
  1. run_fep.py    — pmx hybrid topology + gmx grompp/mdrun for all λ-windows
  2. parse_fep.py  — pymbar MBAR to compute ΔG_apo, ΔG_holo, ΔΔG_bind

PROLINE GATE:
  pmx does not support proline mutations. Sites where ancestral_aa='P' or
  derived_aa='P' are silently skipped by run_fep.py (--skip-proline, default on).
  Skipped sites are written to results/fep/{family}/skipped.tsv.

FEP SCOPE:
  Only run on variants where FoldX × EVmutation quadrant classification is
  "functional_gain" or "stability_function_tradeoff" AND ESM-2 LLR > 0.
  This is enforced through the filtered_variants_tsv input.
  Maximum targets per family: config.fep.max_targets.

Citations:
  Gapsys V, Michielssens S, Seeliger D, de Groot BL. (2015) JCTC 11:4494. [pmx]
  Shirts MR, Chodera JD. (2008) JCP 129:124105. [MBAR]
  Wang L et al. (2015) JACS 137:2695. [RBFE validation benchmark]
"""


rule setup_and_run_fep:
    """Run alchemical FEP for all priority convergent variants in one family.

    Takes the docking-filtered apo and holo PDBs and the quadrant-filtered
    mutation list. Produces one directory tree per variant.
    """
    input:
        apo_pdb="results/foldx/{family}/apo.repaired.pdb",
        holo_pdb="results/docking/{family}/top_pose.pdb",
        filtered_variants="results/integration/{family}.fep_candidates.tsv",
    output:
        sentinel="results/fep/{family}/.setup_done",
    params:
        n_lambda=config["fep"]["n_lambda_windows"],
        ns_per_window=config["fep"]["ns_per_window"],
        n_replicates=config["fep"]["n_replicates"],
        skip_proline="--skip-proline" if config["fep"]["skip_proline"] else "--no-skip-proline",
    conda:
        "../envs/md.yaml"
    log:
        "logs/fep/{family}.setup.log",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=960,           # 16 h; n_lambda × ns × replicates accumulates fast
        gpu=1,
        slurm_partition="gpu",
    shell:
        """
        python workflow/scripts/run_fep.py \
            --apo-pdb        {input.apo_pdb} \
            --holo-pdb       {input.holo_pdb} \
            --mutation-tsv   {input.filtered_variants} \
            --output-dir     results/fep/{wildcards.family} \
            --n-lambda       {params.n_lambda} \
            --ns-per-window  {params.ns_per_window} \
            --n-replicates   {params.n_replicates} \
            {params.skip_proline} \
            --verbose \
            2> {log}

        touch {output.sentinel}
        """


rule parse_fep_results:
    """Run MBAR on all λ-window outputs to compute ΔΔG_bind per variant."""
    input:
        sentinel="results/fep/{family}/.setup_done",
    output:
        tsv="results/fep/{family}.fep_results.tsv",
    conda:
        "../envs/md.yaml"
    log:
        "logs/fep/{family}.mbar.log",
    threads: 2
    resources:
        mem_mb=16000,
        runtime=60,
    shell:
        """
        python workflow/scripts/parse_fep.py \
            --fep-dir  results/fep/{wildcards.family} \
            --output   {output.tsv} \
            --verbose \
            2> {log}
        """
