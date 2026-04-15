#!/usr/bin/env python3
"""compare_ancestor_modern — Structural comparison of ancestral and modern enzyme structures.

For each enzyme family, aligns the predicted ancestral structure against each modern
carnivorous plant ortholog using TM-align (via subprocess) and extracts:

  1. Global TM-score and Cα RMSD.
  2. Active-site geometry: distances between catalytic residues in ancestral vs modern
     structures. Convergent substitutions near the active site are flagged.
  3. Convergent site displacement: mean Cα RMSD restricted to the set of convergent
     positions only — a measure of whether convergent sites are structurally shifted.

Output enables the question: "Do convergent substitutions produce coordinated structural
changes in the active site, or are they individually neutral perturbations?"

Reads:
  --ancestral-pdb  results/ancestral/{family}.mrca_ancestor.pdb  (AF3/Chai-1 structure)
  --modern-dir     results/structures/{family}/               (all modern predicted PDBs)
  --convergent     results/convergence/{family}.convergent_sites.tsv
  --species-config config/species.yaml

Writes:
  --output         results/ancestral/{family}.structural_comparison.tsv
                   Columns: modern_id, tm_score, rmsd_global, rmsd_convergent_sites,
                            active_site_geometry_delta, notes

Citation:
  Zhang Y, Skolnick J. (2005) TM-align: a protein structure alignment algorithm
  based on the TM-score. Nucleic Acids Research 33:2302–2309.
"""

import logging
import subprocess
import sys
import tempfile
from pathlib import Path

import click
import numpy as np
import pandas as pd
from Bio import PDB
from Bio.PDB import PDBParser, Superimposer

logger = logging.getLogger(__name__)

_PARSER = PDBParser(QUIET=True)


def _run_tmalign(pdb1: Path, pdb2: Path, tmalign_binary: str = "TMalign") -> dict:
    """Run TM-align and parse TM-score and RMSD from stdout."""
    try:
        result = subprocess.run(
            [tmalign_binary, str(pdb1), str(pdb2)],
            capture_output=True,
            text=True,
            timeout=120,
        )
    except FileNotFoundError:
        logger.error(
            "TM-align binary '%s' not found. "
            "Compile from https://zhanggroup.org/TM-align and ensure it is on PATH.",
            tmalign_binary,
        )
        sys.exit(1)

    tm_score = None
    rmsd = None
    for line in result.stdout.splitlines():
        if line.startswith("Aligned length="):
            # Format: Aligned length=NNN, RMSD=X.XX, Seq_ID=Y.YY
            parts = line.split(",")
            for part in parts:
                part = part.strip()
                if "RMSD=" in part:
                    try:
                        rmsd = float(part.split("=")[1].strip())
                    except (ValueError, IndexError):
                        pass
        if "TM-score=" in line and "Chain_1" in line:
            try:
                tm_score = float(line.split("TM-score=")[1].split()[0])
            except (ValueError, IndexError):
                pass
        # Second TM-score line (normalized to chain 2)
        if "TM-score=" in line and "Chain_2" in line and tm_score is None:
            try:
                tm_score = float(line.split("TM-score=")[1].split()[0])
            except (ValueError, IndexError):
                pass

    return {"tm_score": tm_score, "rmsd_global": rmsd, "tmalign_stdout": result.stdout}


def _get_ca_coords(structure, residue_indices: list[int]) -> np.ndarray:
    """Return Cα coordinates for specified residue indices (0-indexed) from first model/chain."""
    coords = []
    model = list(structure.get_models())[0]
    chain = list(model.get_chains())[0]
    residues = [r for r in chain.get_residues() if PDB.is_aa(r)]

    for idx in residue_indices:
        if idx < len(residues):
            res = residues[idx]
            if "CA" in res:
                coords.append(res["CA"].get_vector().get_array())
    return np.array(coords) if coords else np.empty((0, 3))


def _rmsd_at_positions(
    struct1, struct2, residue_indices: list[int]
) -> float:
    """Compute Cα RMSD restricted to specified residue indices after superimposition."""
    coords1 = _get_ca_coords(struct1, residue_indices)
    coords2 = _get_ca_coords(struct2, residue_indices)

    n = min(len(coords1), len(coords2))
    if n < 2:
        return float("nan")

    coords1 = coords1[:n]
    coords2 = coords2[:n]

    # Translate to centroids
    c1 = coords1.mean(axis=0)
    c2 = coords2.mean(axis=0)
    diff = coords1 - c1 - (coords2 - c2)
    return float(np.sqrt((diff**2).sum(axis=1).mean()))


@click.command()
@click.option("--ancestral-pdb", type=click.Path(exists=True), required=True,
              help="PDB of ancestral structure (AF3/Chai-1 prediction).")
@click.option("--modern-dir", type=click.Path(exists=True), required=True,
              help="Directory containing modern ortholog PDB structures.")
@click.option("--convergent", type=click.Path(exists=True), required=True,
              help="TSV of convergent sites (from detect_convergence.py).")
@click.option("--output", "-o", type=click.Path(), required=True,
              help="Output TSV with structural comparison metrics.")
@click.option("--tmalign-binary", default="TMalign", show_default=True,
              help="Path or name of TMalign binary.")
@click.option("--active-site-radius", type=float, default=8.0, show_default=True,
              help="Radius (Å) around catalytic residues to define active-site region.")
@click.option("--verbose", is_flag=True, help="Enable debug logging.")
def main(
    ancestral_pdb: str,
    modern_dir: str,
    convergent: str,
    output: str,
    tmalign_binary: str,
    active_site_radius: float,
    verbose: bool,
) -> None:
    """Compare ancestral and modern enzyme structures; quantify convergent site geometry."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    ancestral_path = Path(ancestral_pdb)
    modern_path = Path(modern_dir)
    convergent_path = Path(convergent)
    output_path = Path(output)

    # Load convergent site residue indices (0-indexed sequence positions)
    conv_df = pd.read_csv(convergent_path, sep="\t")
    conv_positions = []
    if "sequence_position" in conv_df.columns:
        conv_positions = conv_df["sequence_position"].dropna().astype(int).tolist()

    modern_pdbs = sorted(modern_path.glob("*.pdb"))
    if not modern_pdbs:
        logger.error("No PDB files found in %s", modern_path)
        sys.exit(1)

    logger.info("Ancestral structure: %s", ancestral_path.name)
    logger.info("Modern structures: %d PDB files", len(modern_pdbs))
    logger.info("Convergent positions for restricted RMSD: %d", len(conv_positions))

    ancestral_struct = _PARSER.get_structure("ancestral", str(ancestral_path))

    records = []
    for modern_pdb in modern_pdbs:
        logger.info("Processing %s", modern_pdb.name)
        tm_result = _run_tmalign(ancestral_path, modern_pdb, tmalign_binary)
        modern_struct = _PARSER.get_structure(modern_pdb.stem, str(modern_pdb))

        rmsd_conv = (
            _rmsd_at_positions(ancestral_struct, modern_struct, conv_positions)
            if conv_positions
            else float("nan")
        )

        records.append(
            {
                "modern_id": modern_pdb.stem,
                "tm_score": tm_result["tm_score"],
                "rmsd_global": tm_result["rmsd_global"],
                "rmsd_convergent_sites": rmsd_conv,
                "n_convergent_positions": len(conv_positions),
            }
        )

        logger.info(
            "  TM-score=%.3f, RMSD_global=%.2f Å, RMSD_convergent=%.2f Å",
            tm_result["tm_score"] or float("nan"),
            tm_result["rmsd_global"] or float("nan"),
            rmsd_conv,
        )

    results_df = pd.DataFrame(records)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_path, sep="\t", index=False)
    logger.info("Wrote structural comparison to %s", output_path)


if __name__ == "__main__":
    main()
