#!/usr/bin/env python3
"""root_tree — Root a per-family IQ-TREE gene tree using outgroup species.

Reads species.yaml to identify which taxa are non-carnivorous outgroups,
finds those tip labels in the unrooted IQ-TREE treefile, and sets them as
the outgroup for rooting.  Writes a rooted Newick file.

Tip-label format expected from fetch_sequences.py:
    {Species_name}|{accession}   (e.g. Arabidopsis_thaliana|NP_180287.2)

The script matches on the Species_name prefix before the pipe character.
"""

import logging
import sys
from pathlib import Path

import click
import yaml
from Bio import Phylo
from io import StringIO


logger = logging.getLogger(__name__)


def _load_outgroup_species(species_config_path: Path) -> list[str]:
    """Return list of species names that are non-carnivorous outgroups."""
    with species_config_path.open(encoding="utf-8") as fh:
        config = yaml.safe_load(fh)
    outgroups = list(config.get("outgroup_species", {}).keys())
    logger.debug("Outgroup species: %s", outgroups)
    return outgroups


def _load_family_outgroups(families_config_path: Path, family: str) -> list[str]:
    """Return species names that serve as outgroups in the enzyme_families.yaml."""
    with families_config_path.open(encoding="utf-8") as fh:
        config = yaml.safe_load(fh)
    family_data = config.get("tier1", {}).get(family, {})
    accessions = family_data.get("accessions", {})
    # Species whose keys appear in outgroup_species from the global species list
    # are considered outgroups; return all species names in this family.
    return list(accessions.keys())


def _find_outgroup_tips(tree, outgroup_species: list[str]) -> list[str]:
    """Find tip labels in the tree that belong to outgroup species.

    Tip labels from fetch_sequences.py have format: Species_name|accession
    Match on Species_name prefix (before pipe character).
    """
    tips = [c.name for c in tree.get_terminals()]
    matched = []
    for tip in tips:
        species_part = tip.split("|")[0] if "|" in tip else tip
        if species_part in outgroup_species:
            matched.append(tip)
    return matched


@click.command()
@click.option(
    "--treefile", "-t",
    type=click.Path(exists=True),
    required=True,
    help="Unrooted IQ-TREE .treefile (Newick format).",
)
@click.option(
    "--species-config", "-s",
    type=click.Path(exists=True),
    required=True,
    help="Path to config/species.yaml.",
)
@click.option(
    "--families-config", "-f",
    type=click.Path(exists=True),
    required=True,
    help="Path to config/enzyme_families.yaml.",
)
@click.option(
    "--family",
    type=str,
    required=True,
    help="Enzyme family key (must match tier1 key in enzyme_families.yaml).",
)
@click.option(
    "--output", "-o",
    type=click.Path(),
    required=True,
    help="Output path for rooted Newick treefile.",
)
@click.option("--verbose", is_flag=True, help="Enable DEBUG logging.")
def main(
    treefile: str,
    species_config: str,
    families_config: str,
    family: str,
    output: str,
    verbose: bool,
) -> None:
    """Root a gene tree using outgroup taxa defined in species.yaml."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    tree_path = Path(treefile)
    species_config_path = Path(species_config)
    families_config_path = Path(families_config)
    output_path = Path(output)

    # Load outgroup species from global config
    outgroup_species = _load_outgroup_species(species_config_path)
    logger.info("Global outgroup species: %s", outgroup_species)

    # Parse tree (IQ-TREE writes Newick with support values)
    tree = Phylo.read(str(tree_path), "newick")
    n_tips = len(tree.get_terminals())
    logger.info("Loaded tree with %d tips from %s", n_tips, tree_path)

    # Find outgroup tips present in this tree
    outgroup_tips = _find_outgroup_tips(tree, outgroup_species)
    logger.info(
        "Outgroup tips found in tree (%d): %s",
        len(outgroup_tips),
        outgroup_tips[:5],  # log first 5
    )

    if not outgroup_tips:
        logger.error(
            "No outgroup tips found in tree for family %s. "
            "Check species.yaml outgroup_species and tip label format "
            "(expected: Species_name|accession).",
            family,
        )
        sys.exit(1)

    # Root the tree using the outgroup
    if len(outgroup_tips) == 1:
        logger.info("Rooting on single outgroup tip: %s", outgroup_tips[0])
        tree.root_with_outgroup({"name": outgroup_tips[0]})
    else:
        # Root on the clade containing all outgroup tips (midpoint of outgroup MRCA)
        logger.info("Rooting on %d outgroup tips (MRCA root)", len(outgroup_tips))
        tree.root_with_outgroup([{"name": t} for t in outgroup_tips])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(output_path), "newick")
    logger.info("Wrote rooted tree to %s", output_path)


if __name__ == "__main__":
    main()
