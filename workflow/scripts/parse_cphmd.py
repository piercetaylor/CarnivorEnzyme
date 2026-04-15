#!/usr/bin/env python3
"""parse_cphmd — Extract pKa values from GROMACS CpHMD λ-dynamics trajectories.

Reads the lambda-state time series from GROMACS CpHMD replica output and computes
pKa for each titratable residue using the Henderson-Hasselbalch relationship between
λ-population at each pH.

The λ coordinate tracks the protonation state: λ ≈ 0 → protonated, λ ≈ 1 → deprotonated.
The protonated fraction = mean(λ < 0.5) over the production trajectory. The pKa is
the pH at which the protonated fraction = 0.5 (interpolated across the pH series).

This script processes one simulation directory (one pH, one variant). Downstream
integration across pH values is done in compare_foldx_evmutation.py.

Reads:
  --replica-dirs  Comma-separated paths to GROMACS replica output directories
                  Each must contain cphmd.xvg (lambda time series from gmx mdrun)

Writes:
  --output  results/cphmd/{family}/{variant}_pH{ph}.pka.tsv
            Columns: residue, residue_number, chain, ph, protonated_fraction,
                     lambda_mean, lambda_std, n_frames

Citation:
  Gapsys V et al. (2022) JCTC 18:4813. GROMACS CpHMD implementation.
  Baptista AM, Soares CM. (2001) Some theoretical and computational aspects of the
  inclusion of proton isomerism in the proton-protein titration. J Phys Chem B.
"""

import logging
import re
import sys
from pathlib import Path

import click
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _parse_lambda_xvg(xvg_path: Path) -> pd.DataFrame:
    """Parse GROMACS .xvg file containing λ time series for all titratable residues.

    Returns DataFrame with columns: time_ps, lambda_1, lambda_2, ... (one per titratable residue).
    """
    lines = []
    headers = []
    with open(xvg_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("#") or line.startswith("@"):
                headers.append(line)
                continue
            cols = line.split()
            if cols:
                lines.append([float(x) for x in cols])

    if not lines:
        logger.warning("No data rows in %s", xvg_path)
        return pd.DataFrame()

    data = np.array(lines)
    col_names = ["time_ps"] + [f"lambda_{i}" for i in range(1, data.shape[1])]
    df = pd.DataFrame(data, columns=col_names)
    return df


def _extract_titratable_residue_info(replica_dirs: list[Path]) -> dict[int, str]:
    """Recover titratable residue labels from phbuilder-generated annotation files.

    phbuilder writes a 'lambda_coupling.dat' file with lines:
        <lambda_index>  <residue_name>  <residue_number>  <chain>

    If that file is absent, falls back to parsing the .top file for TITRATABLE
    residue entries (PHE, ASP, GLU, HIS, LYS, CYS, TYR with alternate protonation).
    Returns {lambda_index (1-based int): "RES_NNN_CHAIN"} or {} on failure.
    """
    # phbuilder writes this in the simulation parent directory (one level up from replicas)
    search_dirs = set()
    for d in replica_dirs:
        search_dirs.add(d.parent)
        search_dirs.add(d.parent.parent)

    for search_dir in search_dirs:
        # Try phbuilder lambda_coupling.dat annotation file
        coupling_file = search_dir / "lambda_coupling.dat"
        if coupling_file.exists():
            logger.debug("Reading titratable residue labels from %s", coupling_file)
            result = {}
            try:
                with open(coupling_file) as fh:
                    for line in fh:
                        line = line.strip()
                        if not line or line.startswith("#"):
                            continue
                        parts = line.split()
                        if len(parts) >= 3:
                            lam_idx = int(parts[0])
                            res_name = parts[1]
                            res_num = parts[2]
                            chain = parts[3] if len(parts) > 3 else "A"
                            result[lam_idx] = f"{res_name}_{res_num}_{chain}"
                if result:
                    logger.info(
                        "Loaded %d titratable residue labels from %s",
                        len(result), coupling_file,
                    )
                    return result
            except (ValueError, IndexError) as exc:
                logger.warning("Could not parse %s: %s", coupling_file, exc)

        # Fallback: scan GROMACS .top for [ moleculetype ] + titratable residue markers
        top_file = search_dir / "titratable.top"
        if top_file.exists():
            logger.debug("Attempting residue label extraction from %s", top_file)
            titratable_aa = {"ASP", "ASPH", "GLU", "GLUH", "HIS", "HISD", "HISE",
                             "LYS", "LYSH", "CYS", "CYSH", "TYR", "TYRH"}
            result = {}
            lam_idx = 1
            try:
                with open(top_file) as fh:
                    for line in fh:
                        parts = line.strip().split()
                        if len(parts) >= 4 and parts[0].upper() in titratable_aa:
                            res_name = parts[0].upper()
                            res_num = parts[2]
                            chain = "A"
                            result[lam_idx] = f"{res_name}_{res_num}_{chain}"
                            lam_idx += 1
                if result:
                    logger.info(
                        "Inferred %d titratable residue labels from topology",
                        len(result),
                    )
                    return result
            except (ValueError, IndexError) as exc:
                logger.warning("Could not parse topology %s: %s", top_file, exc)

    logger.warning(
        "Could not determine titratable residue labels; output will use integer "
        "lambda_index. To fix: ensure phbuilder writes lambda_coupling.dat in %s.",
        [str(d.parent) for d in replica_dirs[:2]],
    )
    return {}


def _compute_protonated_fraction(lambda_series: pd.Series, threshold: float = 0.5) -> tuple[float, float, float]:
    """Return (protonated_fraction, lambda_mean, lambda_std).

    Protonated fraction = fraction of frames where λ < threshold.
    """
    arr = lambda_series.dropna().values
    if len(arr) == 0:
        return float("nan"), float("nan"), float("nan")
    protonated_fraction = float(np.mean(arr < threshold))
    return protonated_fraction, float(np.mean(arr)), float(np.std(arr))


@click.command()
@click.option("--replica-dirs", required=True,
              help="Comma-separated paths to GROMACS CpHMD replica output directories.")
@click.option("--ph", "ph_value", type=float, required=True,
              help="pH value at which this simulation was run.")
@click.option("--variant-id", required=True,
              help="Identifier for this variant (e.g., 'wt' or 'E45Q').")
@click.option("--output", "-o", type=click.Path(), required=True,
              help="Output TSV with per-residue protonation statistics.")
@click.option("--equilibration-ns", type=float, default=10.0, show_default=True,
              help="Equilibration time to discard before analysis (nanoseconds).")
@click.option("--verbose", is_flag=True, help="Enable debug logging.")
def main(
    replica_dirs: str,
    ph_value: float,
    variant_id: str,
    output: str,
    equilibration_ns: float,
    verbose: bool,
) -> None:
    """Parse GROMACS CpHMD λ-dynamics trajectory and compute protonated fractions per residue."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    dir_paths = [Path(d.strip()) for d in replica_dirs.split(",") if d.strip()]
    output_path = Path(output)
    equil_ps = equilibration_ns * 1000.0

    # Attempt to map lambda_index integers to human-readable residue labels
    residue_labels = _extract_titratable_residue_info(dir_paths)

    all_records = []

    for rep_idx, rep_dir in enumerate(dir_paths):
        xvg_path = rep_dir / "cphmd.xvg"
        if not xvg_path.exists():
            logger.warning("Lambda XVG not found: %s — skipping replica %d", xvg_path, rep_idx)
            continue

        df = _parse_lambda_xvg(xvg_path)
        if df.empty:
            logger.warning("Empty lambda data in %s — skipping", xvg_path)
            continue

        # Discard equilibration frames
        df = df[df["time_ps"] >= equil_ps].copy()
        logger.info(
            "Replica %d: %d production frames after discarding %.0f ps equilibration",
            rep_idx, len(df), equil_ps,
        )

        lambda_cols = [c for c in df.columns if c.startswith("lambda_")]
        for col in lambda_cols:
            lambda_idx = int(col.split("_")[1])
            prot_frac, lam_mean, lam_std = _compute_protonated_fraction(df[col])
            all_records.append(
                {
                    "variant_id": variant_id,
                    "ph": ph_value,
                    "replica": rep_idx,
                    "lambda_index": lambda_idx,
                    # residue_label is the human-readable ID (e.g. "GLU_188_A").
                    # Falls back to "lambda_{idx}" if phbuilder annotation is unavailable.
                    "residue_label": residue_labels.get(lambda_idx, f"lambda_{lambda_idx}"),
                    "protonated_fraction": prot_frac,
                    "lambda_mean": lam_mean,
                    "lambda_std": lam_std,
                    "n_frames": len(df),
                }
            )

    if not all_records:
        logger.error("No valid replica data found in any of: %s", [str(d) for d in dir_paths])
        sys.exit(1)

    results_df = pd.DataFrame(all_records)

    # Average across replicas per residue
    agg_df = (
        results_df.groupby(["variant_id", "ph", "lambda_index", "residue_label"])
        .agg(
            protonated_fraction_mean=("protonated_fraction", "mean"),
            protonated_fraction_std=("protonated_fraction", "std"),
            lambda_mean=("lambda_mean", "mean"),
            n_frames_total=("n_frames", "sum"),
            n_replicas=("replica", "count"),
        )
        .reset_index()
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    agg_df.to_csv(output_path, sep="\t", index=False)
    logger.info(
        "Wrote pKa statistics for %d titratable sites to %s",
        agg_df["lambda_index"].nunique(), output_path,
    )


if __name__ == "__main__":
    main()
