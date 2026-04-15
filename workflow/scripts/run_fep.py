#!/usr/bin/env python3
"""run_fep — Alchemical free-energy perturbation for substrate binding ΔΔG.

Computes the change in substrate binding free energy (ΔΔG_bind) caused by each
convergent amino acid substitution using GROMACS + pmx hybrid topology method.

Protocol (pmx alchemical FEP, Gapsys et al. 2015 JCTC):
  1. pmx mutate  → build hybrid topology for ancestral → derived substitution
  2. pmx gentop  → generate GROMACS topology with dual-topology Hamiltonians
  3. gmx grompp  → assemble .tpr for each λ-window (0 = ancestral, 1 = derived)
  4. gmx mdrun   → run short MD at each λ-window; 13 windows × 5 ns × 3 replicates
  5. parse_fep.py → MBAR analysis to compute ΔG and propagate to ΔΔG_bind

Two legs of the thermodynamic cycle:
  - Leg A: mutation in the apo (unbound) protein
  - Leg B: mutation in the holo (substrate-bound) protein
  ΔΔG_bind = ΔG_B − ΔG_A

PROLINE LIMITATION:
  pmx does not support mutations involving proline residues (Pro in either ancestral or
  derived state). Affected convergent sites are skipped with a WARNING log line and
  recorded in a separate skip log. This is not a pipeline error; set config.fep.skip_proline
  to true (the default) to enable automatic skipping.

Reads:
  --apo-pdb       Repaired apo structure (protein only, no substrate)
  --holo-pdb      Repaired holo structure (protein + docked substrate pose)
  --mutation-tsv  TSV with convergent sites: position, ancestral_aa, derived_aa
  --output-dir    Directory for GROMACS FEP run files

Writes:
  --output-dir/{variant}/leg_a/   13 λ-window directories for apo mutation
  --output-dir/{variant}/leg_b/   13 λ-window directories for holo mutation
  --output-dir/skipped.tsv        Sites skipped (proline or other reason)

Citation:
  Gapsys V, Michielssens S, Seeliger D, de Groot BL. (2015) pmx: Automated protein
  structure and topology generation for alchemical perturbations.
  JCTC 11:4494. doi:10.1021/acs.jctc.5b00196
"""

import logging
import subprocess
import sys
from pathlib import Path

import click
import pandas as pd

logger = logging.getLogger(__name__)

# Residues pmx cannot handle in alchemical perturbations
_PMX_UNSUPPORTED = {"PRO", "P"}


def _is_proline_mutation(ancestral_aa: str, derived_aa: str) -> bool:
    """Return True if either amino acid is proline (pmx unsupported)."""
    return ancestral_aa.upper() in _PMX_UNSUPPORTED or derived_aa.upper() in _PMX_UNSUPPORTED


def _check_tools() -> None:
    for binary, name in [("gmx", "GROMACS"), ("pmx", "pmx")]:
        result = subprocess.run(["which", binary], capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(
                "%s ('%s') not found on PATH. "
                "Install: conda install -c conda-forge gromacs; pip install pmx",
                name, binary,
            )
            sys.exit(1)


def _pmx_mutate(pdb_path: Path, position: int, ancestral_aa: str, derived_aa: str, output_dir: Path) -> Path:
    """Run pmx mutate to introduce a hybrid residue at the specified position.

    Returns path to the mutated PDB.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    mutant_pdb = output_dir / f"mutant_{ancestral_aa}{position}{derived_aa}.pdb"

    cmd = [
        "pmx", "mutate",
        "-f", str(pdb_path),
        "-o", str(mutant_pdb),
        "--resnum", str(position),
        "--mut", derived_aa,
        "--ff", "amber99sb-ildn",
    ]
    logger.info("Running pmx mutate: %s%d%s", ancestral_aa, position, derived_aa)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("pmx mutate failed:\n%s", result.stderr[-1500:])
        raise RuntimeError(f"pmx mutate failed for {ancestral_aa}{position}{derived_aa}")

    return mutant_pdb


def _write_lambda_mdp(output_dir: Path, lambda_val: float, window_idx: int, ns_per_window: float) -> Path:
    """Write a GROMACS MDP file for one λ-window of the FEP run."""
    mdp_path = output_dir / f"lambda_{window_idx:02d}.mdp"
    nsteps = int(ns_per_window * 1_000_000 / 2)   # 2 fs timestep

    content = f"""; FEP window {window_idx}: lambda = {lambda_val:.4f}
integrator        = sd
nsteps            = {nsteps}
dt                = 0.002
nstlog            = 5000
nstenergy         = 500
nstxout-compressed = 5000
; Temperature
tc-grps           = System
tau-t             = 2.0
ref-t             = 300
; Pressure
pcoupl            = Parrinello-Rahman
ref-p             = 1.0
compressibility   = 4.5e-5
; PME electrostatics
coulombtype       = PME
rcoulomb          = 1.0
rvdw              = 1.0
; FEP settings
free-energy       = yes
init-lambda-state = {window_idx}
; Lambda vectors (13 windows)
vdw-lambdas       = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.0 1.0
coul-lambdas      = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.5 1.0
mass-lambdas      = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
; dH/dλ output for MBAR
calc-lambda-neighbors = -1
dhdl-print-energy = total
separate-dhdl-file = yes
"""
    mdp_path.write_text(content)
    return mdp_path


def _run_fep_leg(
    pdb_path: Path,
    output_dir: Path,
    n_lambda: int,
    ns_per_window: float,
    n_replicates: int,
    leg_name: str,
    dry_run: bool = False,
) -> None:
    """Set up and run all λ-windows for one leg of the FEP cycle."""
    leg_dir = output_dir / leg_name
    leg_dir.mkdir(parents=True, exist_ok=True)

    lambda_values = [i / (n_lambda - 1) for i in range(n_lambda)]

    for rep in range(n_replicates):
        for win_idx, lam in enumerate(lambda_values):
            win_dir = leg_dir / f"rep_{rep}" / f"window_{win_idx:02d}"
            win_dir.mkdir(parents=True, exist_ok=True)

            mdp_path = _write_lambda_mdp(win_dir, lam, win_idx, ns_per_window)
            tpr_path = win_dir / "run.tpr"

            grompp_cmd = [
                "gmx", "grompp",
                "-f", str(mdp_path),
                "-c", str(pdb_path),
                "-o", str(tpr_path),
                "-maxwarn", "5",
            ]
            logger.debug("grompp: window %d/%d rep %d", win_idx + 1, n_lambda, rep)
            if not dry_run:
                r = subprocess.run(grompp_cmd, capture_output=True, text=True, cwd=str(win_dir))
                if r.returncode != 0:
                    logger.error("grompp failed at window %d rep %d:\n%s", win_idx, rep, r.stderr[-1000:])
                    raise RuntimeError(f"grompp failed: {leg_name} window {win_idx} rep {rep}")

            mdrun_cmd = [
                "gmx", "mdrun",
                "-s", str(tpr_path),
                "-deffnm", "fep",
                "-ntmpi", "1",
                "-ntomp", "4",
            ]
            logger.debug("mdrun: window %d/%d rep %d", win_idx + 1, n_lambda, rep)
            if not dry_run:
                r = subprocess.run(mdrun_cmd, capture_output=True, text=True, cwd=str(win_dir))
                if r.returncode != 0:
                    logger.error("mdrun failed at window %d rep %d:\n%s", win_idx, rep, r.stderr[-1000:])
                    raise RuntimeError(f"mdrun failed: {leg_name} window {win_idx} rep {rep}")

    logger.info(
        "FEP leg '%s' complete: %d windows × %d replicates (dry_run=%s)",
        leg_name, n_lambda, n_replicates, dry_run,
    )


@click.command()
@click.option("--apo-pdb", type=click.Path(exists=True), required=True,
              help="FoldX-repaired apo structure (no substrate).")
@click.option("--holo-pdb", type=click.Path(exists=True), required=True,
              help="FoldX-repaired holo structure (substrate docked).")
@click.option("--mutation-tsv", type=click.Path(exists=True), required=True,
              help="TSV with convergent substitutions to simulate.")
@click.option("--output-dir", "-o", type=click.Path(), required=True,
              help="Root directory for FEP simulation output.")
@click.option("--n-lambda", type=int, default=13, show_default=True,
              help="Number of λ-windows (default 13: standard for MBAR convergence).")
@click.option("--ns-per-window", type=float, default=5.0, show_default=True,
              help="MD time per λ-window in nanoseconds.")
@click.option("--n-replicates", type=int, default=3, show_default=True,
              help="Independent replicates per window for uncertainty estimation.")
@click.option("--skip-proline/--no-skip-proline", default=True, show_default=True,
              help="Skip mutations involving proline (pmx limitation).")
@click.option("--dry-run", is_flag=True,
              help="Write input files only; do not launch GROMACS simulations.")
@click.option("--verbose", is_flag=True, help="Enable debug logging.")
def main(
    apo_pdb: str,
    holo_pdb: str,
    mutation_tsv: str,
    output_dir: str,
    n_lambda: int,
    ns_per_window: float,
    n_replicates: int,
    skip_proline: bool,
    dry_run: bool,
    verbose: bool,
) -> None:
    """Run alchemical FEP for substrate binding ΔΔG at convergent substitution sites."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    apo_path = Path(apo_pdb)
    holo_path = Path(holo_pdb)
    mut_path = Path(mutation_tsv)
    out_path = Path(output_dir)

    if not dry_run:
        _check_tools()

    mutations = pd.read_csv(mut_path, sep="\t")
    required = {"position", "ancestral_aa", "derived_aa"}
    missing = required - set(mutations.columns)
    if missing:
        logger.error("Missing columns in mutation TSV: %s", missing)
        sys.exit(1)

    skipped = []
    ran = []

    for _, row in mutations.iterrows():
        pos = int(row["position"])
        anc = str(row["ancestral_aa"]).upper()
        der = str(row["derived_aa"]).upper()
        variant_id = f"{anc}{pos}{der}"

        if skip_proline and _is_proline_mutation(anc, der):
            reason = "proline_not_supported_by_pmx"
            logger.warning("Skipping %s: %s", variant_id, reason)
            skipped.append({"variant_id": variant_id, "reason": reason})
            continue

        variant_dir = out_path / variant_id
        logger.info("FEP setup for variant %s", variant_id)

        try:
            # Apo leg: mutation in unbound protein
            mutant_apo = _pmx_mutate(apo_path, pos, anc, der, variant_dir / "mutant_prep")
            _run_fep_leg(mutant_apo, variant_dir, n_lambda, ns_per_window, n_replicates, "leg_a", dry_run)

            # Holo leg: mutation in substrate-bound protein
            mutant_holo = _pmx_mutate(holo_path, pos, anc, der, variant_dir / "mutant_prep_holo")
            _run_fep_leg(mutant_holo, variant_dir, n_lambda, ns_per_window, n_replicates, "leg_b", dry_run)

            ran.append({"variant_id": variant_id, "status": "success"})
        except RuntimeError as exc:
            logger.error("FEP failed for %s: %s", variant_id, exc)
            skipped.append({"variant_id": variant_id, "reason": str(exc)})

    # Write skip log
    if skipped:
        skip_path = out_path / "skipped.tsv"
        skip_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(skipped).to_csv(skip_path, sep="\t", index=False)
        logger.info("Skipped %d variants; see %s", len(skipped), skip_path)

    logger.info("FEP complete: %d variants run, %d skipped", len(ran), len(skipped))
    if dry_run:
        logger.info("[DRY RUN] Simulations not launched; input files written.")


if __name__ == "__main__":
    main()
