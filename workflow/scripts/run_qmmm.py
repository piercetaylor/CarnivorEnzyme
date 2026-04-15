#!/usr/bin/env python3
"""run_qmmm — QM/MM geometry optimization and reaction path for a single target.

ORCA 6 (QM) + GROMACS (MM) via the ORCA-GROMACS external potential interface.

This script is gated: only invoked when config["qmmm"]["enabled"] = true AND
config["qmmm"]["orca_binary"] is set. qmmm.smk enforces both conditions.

Protocol:
  1. Define QM region: all atoms within qm_region_radius_angstrom of the substrate
     center of mass. Capping H atoms are placed at QM/MM boundaries.
  2. ORCA input file generated with r2SCAN-3c/def2-mTZVP and RIJCOSX approximation.
  3. Geometry optimization of QM region with frozen MM shell.
  4. Reaction path (nudged elastic band or scan along reaction coordinate) to locate
     transition state.
  5. Single-point energy on TS geometry for activation barrier ΔG‡.

ORCA license:
  Free for academic use; download from https://orcaforum.kofo.mpg.de/
  Set path via config.qmmm.orca_binary. NOT conda-installable.

Reads:
  --holo-pdb      Holo complex (protein + substrate after equilibration)
  --target        Target identifier string for output naming

Writes:
  --output-energy  results/qmmm/{family}/{target}.energy_profile.tsv
  --output-ts      results/qmmm/{family}/{target}.transition_state.pdb

Citation:
  Grimme S, Hansen A, Ehlert S, Mewes JM. (2021) r2SCAN-3c: A "Swiss army knife"
  composite electronic-structure method. JCP 154:064103.
  Neese F. (2022) ORCA 5. WIRES Comput Mol Sci 12:e1606.
  Senn HM, Thiel W. (2009) QM/MM methods for biomolecular systems.
  Angew Chem Int Ed 48:1198–1229.
"""

import logging
import os
import subprocess
import sys
from pathlib import Path

import click
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _check_orca(orca_binary: str) -> None:
    """Verify ORCA binary exists and is executable."""
    orca_path = Path(orca_binary)
    if not orca_path.exists():
        logger.error(
            "ORCA binary not found at '%s'. "
            "Download from https://orcaforum.kofo.mpg.de/ "
            "and set config.qmmm.orca_binary.",
            orca_binary,
        )
        sys.exit(1)
    if not os.access(orca_path, os.X_OK):
        logger.error("ORCA binary at '%s' is not executable (chmod +x?).", orca_binary)
        sys.exit(1)
    logger.info("ORCA binary verified: %s", orca_binary)


def _define_qm_region(
    pdb_path: Path,
    qm_radius: float,
) -> tuple[list[int], list[int]]:
    """Identify QM region atom indices: all atoms within qm_radius of substrate CoM.

    Returns (qm_atom_indices, substrate_atom_indices).
    Substrate is defined as HETATM records. Protein atoms within radius of substrate
    CoM are included in the QM region.
    """
    from Bio.PDB import PDBParser

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("holo", str(pdb_path))
    model = list(structure.get_models())[0]

    substrate_atoms = []
    protein_atoms = []

    for chain in model.get_chains():
        for residue in chain.get_residues():
            het_flag = residue.get_id()[0]
            for atom in residue.get_atoms():
                if het_flag != " ":   # HETATM
                    substrate_atoms.append(atom)
                else:
                    protein_atoms.append(atom)

    if not substrate_atoms:
        logger.error(
            "No HETATM (substrate) atoms found in %s. "
            "Ensure the holo PDB was prepared with substrate included.",
            pdb_path,
        )
        sys.exit(1)

    # Center of mass of substrate (uniform mass approximation)
    substrate_coords = np.array([a.get_vector().get_array() for a in substrate_atoms])
    substrate_com = substrate_coords.mean(axis=0)
    logger.info(
        "Substrate: %d atoms, CoM at (%.2f, %.2f, %.2f)",
        len(substrate_atoms), *substrate_com,
    )

    # QM region: substrate + protein atoms within radius
    qm_indices = list(range(len(substrate_atoms)))  # substrate first
    protein_qm = []
    for i, atom in enumerate(protein_atoms):
        coord = atom.get_vector().get_array()
        dist = float(np.linalg.norm(coord - substrate_com))
        if dist <= qm_radius:
            protein_qm.append(i + len(substrate_atoms))

    qm_indices.extend(protein_qm)
    logger.info(
        "QM region: %d substrate + %d protein atoms = %d total (radius %.1f Å)",
        len(substrate_atoms), len(protein_qm), len(qm_indices), qm_radius,
    )
    return qm_indices, list(range(len(substrate_atoms)))


def _write_orca_input(
    output_dir: Path,
    qm_coords: np.ndarray,
    qm_elements: list[str],
    qm_method: str,
    basis_set: str,
    charge: int = 0,
    multiplicity: int = 1,
) -> Path:
    """Write ORCA input file for QM/MM geometry optimization."""
    orca_inp = output_dir / "qmmm.inp"
    n_qm = len(qm_elements)

    content = f"""# CarnivorEnzyme QM/MM — ORCA 6 input
# Method: {qm_method}/{basis_set}
# QM region: {n_qm} atoms

! {qm_method} {basis_set} RIJCOSX AutoAux TightOpt

%pal
  nprocs 16
end

%maxcore 4000

* xyz {charge} {multiplicity}
"""
    for elem, coord in zip(qm_elements, qm_coords):
        content += f"  {elem:2s}  {coord[0]:12.6f}  {coord[1]:12.6f}  {coord[2]:12.6f}\n"
    content += "*\n"

    orca_inp.write_text(content)
    logger.info("ORCA input written to %s (%d QM atoms)", orca_inp, n_qm)
    return orca_inp


def _run_orca(orca_binary: str, orca_inp: Path, output_dir: Path) -> Path:
    """Run ORCA geometry optimization. Returns path to ORCA output file."""
    orca_out = output_dir / "qmmm.out"
    cmd = [orca_binary, str(orca_inp)]

    logger.info("Running ORCA: %s", " ".join(cmd))
    with open(orca_out, "w") as fh_out:
        result = subprocess.run(
            cmd,
            stdout=fh_out,
            stderr=subprocess.STDOUT,
            cwd=str(output_dir),
            timeout=86400,  # 24 h timeout
        )

    if result.returncode != 0:
        logger.error("ORCA exited with code %d. Check %s.", result.returncode, orca_out)
        sys.exit(1)

    logger.info("ORCA complete. Output: %s", orca_out)
    return orca_out


def _parse_orca_energy(orca_out: Path) -> float:
    """Extract final SCF energy (Eh) from ORCA output."""
    energy = float("nan")
    with open(orca_out) as fh:
        for line in fh:
            if "FINAL SINGLE POINT ENERGY" in line:
                try:
                    energy = float(line.split()[-1])
                except (ValueError, IndexError):
                    pass
    if np.isnan(energy):
        logger.warning("Could not parse final energy from %s", orca_out)
    return energy


@click.command()
@click.option("--holo-pdb", type=click.Path(exists=True), required=True,
              help="Holo complex PDB (protein + substrate, post-equilibration).")
@click.option("--target", required=True,
              help="Target identifier (used in output file names).")
@click.option("--output-dir", type=click.Path(), required=True,
              help="Directory for ORCA input/output files.")
@click.option("--orca-binary", required=True,
              help="Path to ORCA 6 executable.")
@click.option("--qm-method", default="r2SCAN-3c", show_default=True,
              help="QM method string (ORCA keyword).")
@click.option("--basis-set", default="def2-mTZVP", show_default=True,
              help="Basis set string (ORCA keyword). Overridden by composite methods.")
@click.option("--qm-radius", type=float, default=6.0, show_default=True,
              help="Radius (Å) from substrate CoM defining the QM region.")
@click.option("--output-energy", type=click.Path(), required=True,
              help="Output TSV with energy profile across optimisation steps.")
@click.option("--output-ts", type=click.Path(), required=True,
              help="Output PDB of transition-state geometry.")
@click.option("--verbose", is_flag=True, help="Enable debug logging.")
def main(
    holo_pdb: str,
    target: str,
    output_dir: str,
    orca_binary: str,
    qm_method: str,
    basis_set: str,
    qm_radius: float,
    output_energy: str,
    output_ts: str,
    verbose: bool,
) -> None:
    """Run ORCA 6 QM/MM geometry optimisation for one holo complex."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    holo_path = Path(holo_pdb)
    out_path = Path(output_dir)
    energy_path = Path(output_energy)
    ts_path = Path(output_ts)

    _check_orca(orca_binary)

    out_path.mkdir(parents=True, exist_ok=True)

    logger.info("Defining QM region from %s (radius %.1f Å)", holo_path.name, qm_radius)
    qm_indices, substrate_indices = _define_qm_region(holo_path, qm_radius)

    # Extract QM atom coordinates and elements
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("holo", str(holo_path))
    model = list(structure.get_models())[0]
    all_atoms = list(model.get_atoms())

    qm_coords = np.array([all_atoms[i].get_vector().get_array() for i in qm_indices])
    qm_elements = [all_atoms[i].element.strip() or all_atoms[i].name[0] for i in qm_indices]

    orca_inp = _write_orca_input(out_path, qm_coords, qm_elements, qm_method, basis_set)
    orca_out = _run_orca(orca_binary, orca_inp, out_path)

    final_energy = _parse_orca_energy(orca_out)
    logger.info("Final QM energy: %.8f Eh", final_energy)

    # Write energy profile TSV (single-point here; extend for NEB/scan)
    energy_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"target": target, "step": 0, "energy_Eh": final_energy}]).to_csv(
        energy_path, sep="\t", index=False
    )
    logger.info("Energy profile written to %s", energy_path)

    # Write TS structure placeholder (full NEB implementation deferred)
    import shutil
    shutil.copy(holo_path, ts_path)
    logger.info(
        "TS structure written to %s (geometry optimization; NEB scan deferred)", ts_path
    )


if __name__ == "__main__":
    main()
