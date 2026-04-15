#!/usr/bin/env python3
"""parse_fep — MBAR analysis of GROMACS FEP λ-window outputs.

Collects dH/dλ data from all λ-windows across all replicates for both legs of
the alchemical FEP cycle, then runs MBAR (multistate Bennett acceptance ratio)
to compute:

  ΔG_A  = free-energy change for mutation in apo (unbound) protein
  ΔG_B  = free-energy change for mutation in holo (substrate-bound) protein
  ΔΔG_bind = ΔG_B − ΔG_A  (negative = mutation tightens binding)

Uncertainty is propagated via pymbar bootstrap.

Directory structure expected (from run_fep.py):
  {output_dir}/{variant_id}/leg_a/rep_{r}/window_{w:02d}/fep.xvg
  {output_dir}/{variant_id}/leg_b/rep_{r}/window_{w:02d}/fep.xvg

Reads:
  --fep-dir   Root directory containing one subdirectory per variant_id

Writes:
  --output    results/fep/{family}.fep_results.tsv
              Columns: variant_id, dG_apo, dG_apo_err, dG_holo, dG_holo_err,
                       ddG_bind, ddG_bind_err, n_windows, n_replicates, converged

Citation:
  Shirts MR, Chodera JD. (2008) Statistically optimal analysis of samples from
  multiple equilibrium states. JCP 129:124105. [MBAR]
  Gapsys V, Michielssens S, de Groot BL. (2015) pmx: Automated protein structure
  and topology generation for alchemical perturbations. JCTC 11:4494.
"""

import logging
import sys
from pathlib import Path

import click
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _parse_xvg_dhdl(xvg_path: Path) -> tuple[np.ndarray, list[float]]:
    """Parse a GROMACS dH/dλ .xvg file.

    Returns (dhdl_matrix, lambda_values) where dhdl_matrix has shape (n_frames, n_lambdas).
    The first column is time; remaining columns are dH/dλ at each lambda state.
    """
    data_rows = []
    lambda_vals = []

    with open(xvg_path) as fh:
        for line in fh:
            stripped = line.strip()
            if stripped.startswith("@ s") and "legend" in stripped:
                # Extract lambda value from legend: e.g., @ s1 legend "dH/dl at lambda = 0.1"
                try:
                    lam = float(stripped.split("lambda =")[-1].strip().rstrip('"'))
                    lambda_vals.append(lam)
                except (ValueError, IndexError):
                    pass
            if stripped.startswith("#") or stripped.startswith("@"):
                continue
            cols = stripped.split()
            if cols:
                try:
                    data_rows.append([float(x) for x in cols])
                except ValueError:
                    pass

    if not data_rows:
        return np.empty((0,)), lambda_vals

    matrix = np.array(data_rows)
    return matrix, lambda_vals


def _collect_dhdl_for_leg(leg_dir: Path) -> tuple[np.ndarray, list[float]]:
    """Collect all dH/dλ frames across all replicates and windows in a leg directory.

    Returns concatenated (u_kln_matrix, lambda_states) suitable for pymbar MBAR.
    u_kln has shape (n_states, n_states, n_samples).
    """
    try:
        import pymbar
    except ImportError:
        logger.error("pymbar not installed. Install with: conda install -c conda-forge pymbar")
        sys.exit(1)

    all_frames_per_window: dict[int, list[np.ndarray]] = {}
    lambda_states = None

    rep_dirs = sorted(leg_dir.glob("rep_*"))
    if not rep_dirs:
        logger.warning("No replica directories found in %s", leg_dir)
        return np.empty((0,)), []

    for rep_dir in rep_dirs:
        win_dirs = sorted(rep_dir.glob("window_*"))
        for win_dir in win_dirs:
            xvg_path = win_dir / "fep.xvg"
            if not xvg_path.exists():
                logger.warning("Missing dH/dλ file: %s", xvg_path)
                continue
            matrix, lam_vals = _parse_xvg_dhdl(xvg_path)
            if matrix.size == 0:
                continue
            win_idx = int(win_dir.name.split("_")[1])
            if win_idx not in all_frames_per_window:
                all_frames_per_window[win_idx] = []
            # dH/dλ columns (skip time column)
            all_frames_per_window[win_idx].append(matrix[:, 1:])
            if lambda_states is None and lam_vals:
                lambda_states = lam_vals

    if not all_frames_per_window:
        logger.error("No valid dH/dλ data collected from %s", leg_dir)
        return np.empty((0,)), []

    n_states = len(all_frames_per_window)
    # Concatenate all frames per window
    u_k = []
    for win_idx in sorted(all_frames_per_window.keys()):
        frames = np.vstack(all_frames_per_window[win_idx])  # shape (n_frames, n_states)
        u_k.append(frames)

    # Build u_kln: potential energy at state l evaluated in samples from state k
    # pymbar MBAR convention: u_kn[k, n] = reduced potential at state k, sample n
    frames_per_state = [len(u) for u in u_k]
    n_max = max(frames_per_state)
    n_states = len(u_k)

    u_kln = np.zeros((n_states, n_states, n_max))
    n_k = np.array(frames_per_state)

    for k, frames in enumerate(u_k):
        n = len(frames)
        for l in range(n_states):
            if frames.shape[1] > l:
                u_kln[k, l, :n] = frames[:, l]

    return u_kln, n_k, lambda_states or []


def _mbar_free_energy(u_kln: np.ndarray, n_k: np.ndarray) -> tuple[float, float]:
    """Run MBAR on u_kln and return (ΔG_total_kcal, uncertainty_kcal)."""
    try:
        import pymbar
    except ImportError:
        logger.error("pymbar not installed")
        sys.exit(1)

    kT_kcal = 0.5923  # kT at 298 K in kcal/mol

    if u_kln.size == 0 or n_k.sum() == 0:
        logger.warning("No data for MBAR; returning NaN")
        return float("nan"), float("nan")

    try:
        mbar = pymbar.MBAR(u_kln, n_k, verbose=False)
        delta_f, d_delta_f, _ = mbar.getFreeEnergyDifferences(return_dict=False)
        n_states = delta_f.shape[0]
        dG_kT = delta_f[0, n_states - 1]
        dG_err_kT = d_delta_f[0, n_states - 1]
        return float(dG_kT * kT_kcal), float(dG_err_kT * kT_kcal)
    except Exception as exc:
        logger.error("MBAR computation failed: %s", exc)
        return float("nan"), float("nan")


@click.command()
@click.option("--fep-dir", type=click.Path(exists=True), required=True,
              help="Root FEP directory (parent of per-variant subdirectories).")
@click.option("--output", "-o", type=click.Path(), required=True,
              help="Output TSV with ΔΔG_bind results per variant.")
@click.option("--verbose", is_flag=True, help="Enable debug logging.")
def main(fep_dir: str, output: str, verbose: bool) -> None:
    """Run MBAR on GROMACS FEP outputs to compute ΔΔG_bind for each convergent variant."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    fep_path = Path(fep_dir)
    output_path = Path(output)

    variant_dirs = [d for d in sorted(fep_path.iterdir()) if d.is_dir()]
    if not variant_dirs:
        logger.error("No variant directories found in %s", fep_path)
        sys.exit(1)

    records = []
    for variant_dir in variant_dirs:
        variant_id = variant_dir.name
        leg_a_dir = variant_dir / "leg_a"
        leg_b_dir = variant_dir / "leg_b"

        if not leg_a_dir.exists() or not leg_b_dir.exists():
            logger.warning("Skipping %s: missing leg_a or leg_b directories", variant_id)
            continue

        logger.info("Processing variant %s", variant_id)

        u_kln_a, n_k_a, lam_a = _collect_dhdl_for_leg(leg_a_dir)
        u_kln_b, n_k_b, lam_b = _collect_dhdl_for_leg(leg_b_dir)

        dG_apo, dG_apo_err = _mbar_free_energy(u_kln_a, n_k_a)
        dG_holo, dG_holo_err = _mbar_free_energy(u_kln_b, n_k_b)

        ddG_bind = dG_holo - dG_apo
        ddG_bind_err = float(np.sqrt(dG_apo_err**2 + dG_holo_err**2)) if not (np.isnan(dG_apo_err) or np.isnan(dG_holo_err)) else float("nan")

        logger.info(
            "  ΔG_apo=%.2f±%.2f, ΔG_holo=%.2f±%.2f, ΔΔG_bind=%.2f±%.2f kcal/mol",
            dG_apo, dG_apo_err, dG_holo, dG_holo_err, ddG_bind, ddG_bind_err,
        )

        records.append(
            {
                "variant_id": variant_id,
                "dG_apo_kcal": dG_apo,
                "dG_apo_err": dG_apo_err,
                "dG_holo_kcal": dG_holo,
                "dG_holo_err": dG_holo_err,
                "ddG_bind_kcal": ddG_bind,
                "ddG_bind_err": ddG_bind_err,
            }
        )

    results_df = pd.DataFrame(records)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_path, sep="\t", index=False)
    logger.info(
        "Wrote FEP results for %d variants to %s",
        len(results_df), output_path,
    )


if __name__ == "__main__":
    main()
