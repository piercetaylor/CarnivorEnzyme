"""qmmm.smk — Phase 5E: QM/MM for catalytic mechanism of highest-priority targets.

ORCA 6 (QM) + GROMACS (MM) with the ORCA-GROMACS interface. Reserved for a
maximum of 2 targets, selected post-FEP as those with the most compelling ΔΔG_bind
evidence for mechanistic change.

QM/MM is NOT run by default. Two gates must both pass:
  1. config["qmmm"]["enabled"] must be true (set manually after FEP results reviewed)
  2. config["qmmm"]["orca_binary"] must be set to the ORCA 6 executable path

Scientific rationale:
  FEP gives the thermodynamic answer (ΔΔG_bind). QM/MM gives the mechanistic answer:
  what bond is broken/formed, what is the transition state geometry, does the convergent
  substitution lower the activation barrier? This is the level of evidence needed to
  claim a mechanistic (not just affinity) change.

QM method: r2SCAN-3c composite DFT (Grimme et al. 2021, JCP 154:064103)
  - Validated for peptide-enzyme interactions and metalloenzyme active sites
  - ~100× faster than CCSD(T) with acceptable accuracy (MAE ~1.5 kcal/mol vs DLPNO)

ORCA LICENSE:
  ORCA 6 is free for academic use but NOT conda-installable.
  Download: https://orcaforum.kofo.mpg.de/
  Set path in config/config.yaml: qmmm.orca_binary

Citations:
  Grimme S, Hansen A, Ehlert S, Mewes JM. (2021) r2SCAN-3c: A "Swiss army knife"
  composite electronic-structure method. JCP 154:064103.
  Senn HM, Thiel W. (2009) QM/MM methods for biomolecular systems.
  Angew Chem Int Ed 48:1198–1229.
  Neese F. (2022) Software update: The ORCA program system — Version 5.0.
  WIRES Comput Mol Sci 12:e1606.
"""

import os


# ── Gate check ────────────────────────────────────────────────────────────────
# This rule will not be reachable unless qmmm.enabled = true in config.
# Wrap in a conditional so the DAG does not include these rules by default.

if config.get("qmmm", {}).get("enabled", False):

    rule run_qmmm:
        """QM/MM geometry optimization + reaction path for a single target.

        Input PDB must be the holo complex (protein + substrate in binding site)
        after equilibration. The QM region is defined as all atoms within
        config.qmmm.qm_region_radius_angstrom of the substrate.
        """
        input:
            holo_pdb="results/fep/{family}/top_target.pdb",
        output:
            energy_profile="results/qmmm/{family}/{target}.energy_profile.tsv",
            ts_structure="results/qmmm/{family}/{target}.transition_state.pdb",
        params:
            orca_binary=config["qmmm"]["orca_binary"],
            qm_method=config["qmmm"]["qm_method"],
            basis_set=config["qmmm"]["basis_set"],
            qm_radius=config["qmmm"]["qm_region_radius_angstrom"],
        conda:
            "../envs/qmmm.yaml"
        log:
            "logs/qmmm/{family}_{target}.log",
        threads: 16
        resources:
            mem_mb=64000,
            runtime=1440,      # 24 h; QM/MM is expensive
            slurm_partition="general",
        shell:
            """
            # Verify ORCA binary exists (config gate)
            if [ -z "{params.orca_binary}" ] || [ ! -x "{params.orca_binary}" ]; then
                echo "ERROR: ORCA binary not found at '{params.orca_binary}'." >&2
                echo "Download from https://orcaforum.kofo.mpg.de/ and set qmmm.orca_binary." >&2
                exit 1
            fi
            export ORCA_PATH=$(dirname {params.orca_binary})

            python workflow/scripts/run_qmmm.py \
                --holo-pdb       {input.holo_pdb} \
                --target         {wildcards.target} \
                --output-dir     results/qmmm/{wildcards.family}/{wildcards.target} \
                --orca-binary    {params.orca_binary} \
                --qm-method      {params.qm_method} \
                --basis-set      {params.basis_set} \
                --qm-radius      {params.qm_radius} \
                --output-energy  {output.energy_profile} \
                --output-ts      {output.ts_structure} \
                --verbose \
                2> {log}
            """

else:
    # Provide a no-op rule so Snakemake does not fail on include
    rule run_qmmm:
        """QM/MM disabled (config.qmmm.enabled = false). No-op rule."""
        run:
            pass
