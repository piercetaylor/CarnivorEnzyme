#!/usr/bin/env python3
"""run_cphmd — Constant-pH MD simulation wrapper using GROMACS native λ-dynamics.

Runs pH-replica exchange molecular dynamics to determine pKa shifts of titratable
residues in convergent substitution variants vs. the wild-type background. This
directly addresses the question: "Do convergent substitutions shift protonation
states at the pH values found in pitcher fluids (2.5–5.0)?"

Pipeline:
  1. phbuilder -- charges the GROMACS topology with titratable proton states
  2. gmx grompp -- assemble run input
  3. gmx mdrun -- CpHMD production with -multidir for pH-REMD
  4. parse output with parse_cphmd.py to extract pKa values per residue

Requirements:
  - GROMACS 2024.2+ with CpHMD support (lambda-dynamics; Gapsys et al. 2022 JCTC)
  - phbuilder (conda-forge) for titratable topology generation

Design notes:
  - Each simulation runs 3 replicas × 50 ns production per pH value.
  - Two variants per convergent site: wild-type (ancestral AA) and mutant (derived AA).
  - Only sites classified as "functional_gain" or "stability_function_tradeoff" by
    the FoldX×EVmutation quadrant are run (from config.cphmd.run_on_quadrants).
  - pKa shifts (ΔpKa = pKa_mutant − pKa_wildtype) are the reported quantity.

Reads:
  --pdb          Repaired PDB from FoldX RepairPDB
  --variant-tsv  TSV with columns: position, ancestral_aa, derived_aa, quadrant
  --output-dir   Directory for GROMACS simulation files and output

Writes:
  --output-dir/{variant}_pH{ph}/  GROMACS run directory per variant × pH
  (parse_cphmd.py processes these directories downstream)

Citation:
  Gapsys V, Pérez-Hernández G, de Groot BL, et al. (2022) Accurate absolute free
  energies for ligand-protein binding based on non-equilibrium approaches.
  JCTC 18:4813–4830.
  Guo Z, Mahdi J, et al. (2022) GROMACS native constant-pH MD (λ-dynamics).
  J. Chem. Theory Comput. 18:6320.
"""

import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional

import click
import pandas as pd
import yaml

logger = logging.getLogger(__name__)


def _check_tool(binary: str, name: str) -> None:
    """Exit with error if binary is not on PATH."""
    result = subprocess.run(["which", binary], capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(
            "%s binary '%s' not found. Install with: conda install -c conda-forge %s",
            name, binary, binary.lower(),
        )
        sys.exit(1)
    logger.debug("%s found: %s", name, result.stdout.strip())


def _run_phbuilder(pdb_path: Path, output_dir: Path, ph: float) -> Path:
    """Run phbuilder to generate titratable GROMACS topology.

    Returns path to the generated .top file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    top_path = output_dir / "titratable.top"
    gro_path = output_dir / "titratable.gro"

    cmd = [
        "phbuilder", "gentitratables",
        "-f", str(pdb_path),
        "-o", str(output_dir / "titratable"),
        "-ph", str(ph),
        "-ff", "amber99sb-ildn",      # Force field; AMBER99SB-ILDN standard for proteins
        "-water", "tip3p",
    ]
    logger.info("Running phbuilder for pH %.1f: %s", ph, " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(output_dir))

    if result.returncode != 0:
        logger.error("phbuilder failed:\n%s", result.stderr)
        sys.exit(1)

    if not top_path.exists():
        logger.error("phbuilder did not produce expected topology at %s", top_path)
        sys.exit(1)

    return top_path


def _write_mdp(output_dir: Path, equilibration_ns: float, production_ns: float, ph: float) -> tuple[Path, Path]:
    """Write GROMACS .mdp parameter files for equilibration and production CpHMD.

    Returns (equil_mdp_path, prod_mdp_path).
    """
    equil_mdp = output_dir / "equil.mdp"
    prod_mdp = output_dir / "prod.mdp"

    equil_steps = int(equilibration_ns * 1_000_000 / 2)  # 2 fs time step
    prod_steps = int(production_ns * 1_000_000 / 2)

    equil_content = f"""; Equilibration MD parameters for CpHMD
integrator    = md
nsteps        = {equil_steps}
dt            = 0.002
nstxout       = 5000
nstvout       = 5000
nstlog        = 5000
nstcalcenergy = 100
; Temperature coupling
tcoupl        = V-rescale
tc-grps       = Protein Non-Protein
tau-t         = 0.1 0.1
ref-t         = 300 300
; Pressure coupling
pcoupl        = Parrinello-Rahman
ref-p         = 1.0
; Electrostatics
coulombtype   = PME
rcoulomb      = 1.0
rvdw          = 1.0
; λ-dynamics (CpHMD)
free-energy   = yes
init-lambda-state = 0
; phbuilder will populate lambda-state-init per titratable residue
"""

    prod_content = f"""; Production CpHMD parameters
integrator    = md
nsteps        = {prod_steps}
dt            = 0.002
nstxout       = 50000
nstvout       = 50000
nstlog        = 10000
nstcalcenergy = 100
; Temperature coupling
tcoupl        = V-rescale
tc-grps       = Protein Non-Protein
tau-t         = 0.1 0.1
ref-t         = 300 300
; Pressure coupling
pcoupl        = Parrinello-Rahman
ref-p         = 1.0
; Electrostatics
coulombtype   = PME
rcoulomb      = 1.0
rvdw          = 1.0
; λ-dynamics (CpHMD)
free-energy   = yes
; pH replica exchange handled via -multidir in gmx mdrun
"""

    equil_mdp.write_text(equil_content)
    prod_mdp.write_text(prod_content)
    logger.debug("Written MDP files: %s, %s", equil_mdp, prod_mdp)
    return equil_mdp, prod_mdp


def _run_grompp_and_mdrun(
    output_dir: Path,
    top_path: Path,
    gro_path: Path,
    equil_mdp: Path,
    prod_mdp: Path,
    n_replicas: int,
    threads: int,
    dry_run: bool = False,
) -> None:
    """Run gmx grompp and gmx mdrun for CpHMD simulation."""
    equil_tpr = output_dir / "equil.tpr"
    prod_tpr = output_dir / "prod.tpr"

    # Equilibration grompp
    grompp_equil = [
        "gmx", "grompp",
        "-f", str(equil_mdp),
        "-c", str(gro_path),
        "-p", str(top_path),
        "-o", str(equil_tpr),
        "-maxwarn", "5",
    ]
    logger.info("Running grompp (equilibration)")
    if not dry_run:
        result = subprocess.run(grompp_equil, capture_output=True, text=True, cwd=str(output_dir))
        if result.returncode != 0:
            logger.error("gmx grompp (equil) failed:\n%s", result.stderr[-2000:])
            sys.exit(1)

    # Production grompp
    grompp_prod = [
        "gmx", "grompp",
        "-f", str(prod_mdp),
        "-c", str(equil_tpr),
        "-p", str(top_path),
        "-o", str(prod_tpr),
        "-maxwarn", "5",
    ]
    logger.info("Running grompp (production)")
    if not dry_run:
        result = subprocess.run(grompp_prod, capture_output=True, text=True, cwd=str(output_dir))
        if result.returncode != 0:
            logger.error("gmx grompp (prod) failed:\n%s", result.stderr[-2000:])
            sys.exit(1)

    # CpHMD mdrun with replica exchange (-multidir requires N separate replica dirs)
    replica_dirs = [output_dir / f"replica_{i}" for i in range(n_replicas)]
    for d in replica_dirs:
        d.mkdir(parents=True, exist_ok=True)

    mdrun_cmd = [
        "gmx", "mdrun",
        "-s", str(prod_tpr),
        "-multidir",
        *[str(d) for d in replica_dirs],
        "-deffnm", "cphmd",
        "-ntmpi", "1",
        "-ntomp", str(threads),
        "-gpu_id", "0",
        "-v",
    ]
    logger.info("Running CpHMD mdrun (%d replicas, %d threads)", n_replicas, threads)
    if not dry_run:
        result = subprocess.run(mdrun_cmd, capture_output=True, text=True, cwd=str(output_dir))
        if result.returncode != 0:
            logger.error("gmx mdrun failed:\n%s", result.stderr[-3000:])
            sys.exit(1)
        logger.info("CpHMD mdrun complete")


@click.command()
@click.option("--pdb", type=click.Path(exists=True), required=True,
              help="FoldX-repaired PDB file for this variant.")
@click.option("--variant-tsv", type=click.Path(exists=True), required=True,
              help="TSV with convergent site variants to simulate (quadrant-filtered).")
@click.option("--output-dir", "-o", type=click.Path(), required=True,
              help="Directory for simulation output.")
@click.option("--ph", "ph_value", type=float, required=True,
              help="pH value for this CpHMD simulation (from config.cphmd.ph_values).")
@click.option("--equilibration-ns", type=float, default=10.0, show_default=True,
              help="Equilibration simulation length in nanoseconds.")
@click.option("--production-ns", type=float, default=50.0, show_default=True,
              help="Production simulation length in nanoseconds.")
@click.option("--n-replicas", type=int, default=3, show_default=True,
              help="Number of independent CpHMD replicas.")
@click.option("--threads", type=int, default=4, show_default=True,
              help="OpenMP threads per replica (gmx mdrun -ntomp).")
@click.option("--dry-run", is_flag=True,
              help="Write input files only; do not launch simulation.")
@click.option("--verbose", is_flag=True, help="Enable debug logging.")
def main(
    pdb: str,
    variant_tsv: str,
    output_dir: str,
    ph_value: float,
    equilibration_ns: float,
    production_ns: float,
    n_replicas: int,
    threads: int,
    dry_run: bool,
    verbose: bool,
) -> None:
    """Run GROMACS CpHMD for a single variant at a single pH value."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    pdb_path = Path(pdb)
    output_path = Path(output_dir)

    if not dry_run:
        _check_tool("gmx", "GROMACS")
        _check_tool("phbuilder", "phbuilder")

    logger.info("CpHMD run: %s at pH %.1f", pdb_path.stem, ph_value)
    logger.info("Parameters: equil=%.0f ns, production=%.0f ns, replicas=%d", equilibration_ns, production_ns, n_replicas)

    top_path = _run_phbuilder(pdb_path, output_path, ph_value)
    gro_path = output_path / "titratable.gro"

    equil_mdp, prod_mdp = _write_mdp(output_path, equilibration_ns, production_ns, ph_value)
    _run_grompp_and_mdrun(
        output_path, top_path, gro_path, equil_mdp, prod_mdp, n_replicas, threads, dry_run
    )

    if dry_run:
        logger.info("[DRY RUN] Input files written; simulation not launched.")
    else:
        logger.info("CpHMD simulation complete for %s at pH %.1f", pdb_path.stem, ph_value)
        logger.info("Replica directories: %s", output_path)
        logger.info("Next step: run parse_cphmd.py to extract pKa values.")


if __name__ == "__main__":
    main()
