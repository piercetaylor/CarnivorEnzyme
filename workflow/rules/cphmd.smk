"""cphmd.smk — Phase 5C: Constant-pH MD to quantify protonation state shifts.

GROMACS native λ-dynamics CpHMD (Gapsys et al. 2022, JCTC) identifies whether
convergent substitutions shift pKa values of nearby titratable residues. This is
particularly relevant for enzymes active at pitcher fluid pH (2.5–5.0), where
partial protonation of catalytic residues (Asp/Glu) directly modulates activity.

Question: Do convergent substitutions shift the protonation equilibrium of active-site
Asp/Glu residues toward the protonation state optimal for their lineage's pitcher pH?

Pipeline per variant (one convergent substitution × one pH value):
  1. run_cphmd.py      — phbuilder topology setup + GROMACS mdrun (50 ns × 3 replicas)
  2. parse_cphmd.py    — extract λ time series → protonated fraction per residue
  3. Aggregation across pH values → pKa via Henderson-Hasselbalch interpolation
     (this aggregation is done downstream in integrate.smk / build_atlas.py)

Gate: Only sites classified as "functional_gain" or "stability_function_tradeoff"
by the FoldX × EVmutation quadrant (from integrate.smk) are run. This gate is
implemented by the variant_tsv input, which is filtered before this rule fires.

Resource note: 50 ns × 3 replicas × 3 pH values × N variants = substantial compute.
This rule submits to GPU nodes on Hellbender. Expected runtime ~8 h per variant × pH.

Citations:
  Gapsys V, Pérez-Hernández G, de Groot BL, et al. (2022) JCTC 18:4813.
  Guo Z et al. (2022) GROMACS CpHMD. JCTC 18:6320.
  Baptista AM, Soares CM. (2001) J Phys Chem B 105:293.
"""


rule run_cphmd:
    """GROMACS CpHMD simulation for one variant at one pH value.

    Each (family, variant, ph) combination is a separate job. Wildcard ph is
    a string like '2_5' representing pH 2.5 (dot replaced with underscore for
    safe file paths; run_cphmd.py receives the float value via params).
    """
    input:
        pdb="results/foldx/{family}/{variant}.repaired.pdb",
        variant_tsv="results/integration/{family}.quadrant_filtered_variants.tsv",
    output:
        sentinel="results/cphmd/{family}/{variant}/pH{ph}/.done",
    params:
        ph_float=lambda wc: float(wc.ph.replace("_", ".")),
        equilibration_ns=config["cphmd"]["equilibration_ns"],
        production_ns=config["cphmd"]["production_ns"],
        n_replicas=config["cphmd"]["n_replicas"],
    conda:
        "../envs/md.yaml"
    log:
        "logs/cphmd/{family}_{variant}_pH{ph}.log",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=600,          # 10 h; extend to 960 for 100 ns production
        gpu=1,
        slurm_partition="gpu",
    shell:
        """
        python workflow/scripts/run_cphmd.py \
            --pdb             {input.pdb} \
            --variant-tsv     {input.variant_tsv} \
            --output-dir      results/cphmd/{wildcards.family}/{wildcards.variant}/pH{wildcards.ph} \
            --ph              {params.ph_float} \
            --equilibration-ns {params.equilibration_ns} \
            --production-ns   {params.production_ns} \
            --n-replicas      {params.n_replicas} \
            --threads         {threads} \
            --verbose \
            2> {log}

        # Touch sentinel on success
        touch {output.sentinel}
        """


rule parse_cphmd:
    """Parse GROMACS CpHMD λ-trajectories to extract protonated fractions.

    Aggregates across replicas within a single (variant, pH) combination.
    Cross-pH aggregation (pKa fitting) is done in build_atlas.py.
    """
    input:
        sentinel="results/cphmd/{family}/{variant}/pH{ph}/.done",
    output:
        tsv="results/cphmd/{family}/{variant}_pH{ph}.pka.tsv",
    params:
        ph_float=lambda wc: float(wc.ph.replace("_", ".")),
        equilibration_ns=config["cphmd"]["equilibration_ns"],
    conda:
        "../envs/md.yaml"
    log:
        "logs/cphmd/{family}_{variant}_pH{ph}.parse.log",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    shell:
        """
        # Collect replica directories created by run_cphmd.py
        REPLICA_DIRS=$(python -c "
import pathlib
dirs = sorted(pathlib.Path('results/cphmd/{wildcards.family}/{wildcards.variant}/pH{wildcards.ph}').glob('replica_*'))
print(','.join(str(d) for d in dirs))
")

        python workflow/scripts/parse_cphmd.py \
            --replica-dirs    "$REPLICA_DIRS" \
            --ph              {params.ph_float} \
            --variant-id      {wildcards.variant} \
            --output          {output.tsv} \
            --equilibration-ns {params.equilibration_ns} \
            --verbose \
            2> {log}
        """
