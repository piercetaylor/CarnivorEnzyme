#!/usr/bin/env python3
"""extract_ancestor — Parse IQ-TREE .state file and write ancestral FASTA for AF3 prediction.

IQ-TREE with -ancestral writes a .state file containing marginal posterior probabilities
for every ancestral node at every alignment position. This script:

  1. Reads the .state file from IQ-TREE3 (-ancestral output).
  2. Identifies the MRCA (most recent common ancestor) node for a specified set
     of carnivorous plant taxa — the node whose sequence represents the last common
     ancestor BEFORE the carnivorous lineage divergence.
  3. Builds a consensus ancestral sequence by taking the MAP (maximum a posteriori)
     amino acid at each position, subject to a posterior probability threshold.
  4. Writes the ancestral FASTA for downstream structure prediction with AF3/Chai-1.

.state file format (IQ-TREE3):
  Tab-separated, columns: Node, Site, State, {AA probabilities in alphabetical order}
  One row per (node, site) pair. Site is 1-indexed. Header line present.

Reads:
  --state-file    results/phylogenies/{family}.state   (IQ-TREE ancestral .state file)
  --tree-file     results/phylogenies/{family}.rooted.treefile
  --outgroup-ids  Comma-separated list of sequence IDs that are outgroup (non-carnivorous)

Writes:
  --output        results/ancestral/{family}.mrca_ancestor.fa
                  Single FASTA record: ancestral MAP sequence for the carnivore MRCA
  --output-stats  results/ancestral/{family}.mrca_stats.tsv
                  Posterior probability per site; flags uncertain positions

Citation:
  Pupko T, Pe I, Shamir R, Graur D. (2000) A fast algorithm for joint reconstruction
  of ancestral amino acid sequences. Mol Biol Evol 17:890–896.
  Yang Z. (1997) PAML: a program package for phylogenetic analysis by maximum likelihood.
  Comput Appl Biosci 13:555–556.
"""

import logging
import sys
from pathlib import Path
from typing import Optional

import click
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

# Standard amino acid order used by IQ-TREE in .state file columns
_IQTREE_AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

# MRCA node identification uses ete3 tree traversal
_CARNIVOROUS_TAXA_KEYWORDS = [
    "Nepenthes", "Drosera", "Dionaea", "Aldrovanda", "Drosophyllum",
    "Cephalotus", "Sarracenia", "Utricularia", "Pinguicula",
]


def _load_state_file(state_path: Path) -> pd.DataFrame:
    """Load IQ-TREE .state file into a DataFrame.

    Returns DataFrame with columns: Node, Site, State, plus one column per AA.
    """
    logger.info("Reading IQ-TREE .state file from %s", state_path)
    df = pd.read_csv(state_path, sep="\t", comment="#")
    # IQ-TREE .state header: Node  Site  State  p_A  p_C  p_D  ...  p_Y
    # The AA columns may be labeled 'p_A' etc. or just 'A' etc. depending on version.
    logger.info("Loaded %d rows from .state file", len(df))
    return df


def _find_mrca_node(tree_path: Path, taxa_keywords: list[str]) -> str:
    """Use ete3 to find the MRCA node name for all carnivorous taxa.

    Returns the node name string as it appears in the .state file.
    """
    try:
        from ete3 import Tree
    except ImportError:
        logger.error("ete3 not installed. Install with: conda install -c conda-forge ete3")
        sys.exit(1)

    logger.info("Loading tree from %s", tree_path)
    tree = Tree(str(tree_path), format=1)

    # Collect leaf names matching any carnivorous keyword
    carnivorous_leaves = [
        leaf.name for leaf in tree.get_leaves()
        if any(kw.lower() in leaf.name.lower() for kw in taxa_keywords)
    ]

    if not carnivorous_leaves:
        logger.error(
            "No carnivorous taxa found in tree. "
            "Leaf names: %s", [l.name for l in tree.get_leaves()][:10]
        )
        sys.exit(1)

    logger.info("Carnivorous taxa found in tree: %s", carnivorous_leaves)

    if len(carnivorous_leaves) == 1:
        # Only one carnivorous taxon — use the leaf itself as "MRCA"
        mrca = tree.search_nodes(name=carnivorous_leaves[0])[0]
    else:
        mrca = tree.get_common_ancestor(carnivorous_leaves)

    # Get ete3 node name — IQ-TREE labels internal nodes as "NodeX"
    node_name = mrca.name
    if not node_name:
        # If unnamed, IQ-TREE assigns labels like "Node1", "Node2" by traversal order.
        # Use ete3 to identify — fall back to index-based name.
        for i, node in enumerate(tree.traverse("levelorder")):
            if node is mrca:
                node_name = f"Node{i + 1}"
                break

    logger.info("MRCA node identified: %s (spans %d carnivorous taxa)", node_name, len(carnivorous_leaves))
    return node_name


def _build_ancestral_sequence(
    state_df: pd.DataFrame,
    node_name: str,
    posterior_threshold: float,
    aa_col_prefix: str = "p_",
) -> tuple[str, pd.DataFrame]:
    """Build MAP ancestral sequence from .state file rows for a given node.

    Returns (sequence_string, stats_dataframe).
    Stats columns: site, map_aa, map_posterior, uncertain (bool).
    """
    node_df = state_df[state_df["Node"] == node_name].copy()
    if node_df.empty:
        # Try alternate IQ-TREE naming conventions
        alt_name = node_name.replace("Node", "node").replace("node", "Node")
        node_df = state_df[state_df["Node"] == alt_name].copy()
    if node_df.empty:
        # Last resort: partial match (IQ-TREE3 may prefix with '/')
        node_df = state_df[state_df["Node"].str.contains(node_name, na=False)].copy()

    if node_df.empty:
        available = state_df["Node"].unique()[:10].tolist()
        logger.error(
            "Node '%s' not found in .state file. Available nodes (first 10): %s",
            node_name, available,
        )
        sys.exit(1)

    node_df = node_df.sort_values("Site")

    # Detect AA probability column names
    aa_cols = [c for c in node_df.columns if c in _IQTREE_AA_ORDER or c.startswith(aa_col_prefix)]
    if not aa_cols:
        # IQ-TREE3 uses alphabetical single-letter columns directly
        aa_cols = [c for c in node_df.columns if c in _IQTREE_AA_ORDER]

    if not aa_cols:
        logger.error(
            "Cannot identify amino acid probability columns. "
            "Column names found: %s", node_df.columns.tolist()
        )
        sys.exit(1)

    logger.info("Using AA probability columns: %s", aa_cols[:5], "... (%d total)" if len(aa_cols) > 5 else "")

    # Map each column to its single-letter AA
    col_to_aa = {}
    for col in aa_cols:
        stripped = col.replace(aa_col_prefix, "").strip()
        if stripped in _IQTREE_AA_ORDER:
            col_to_aa[col] = stripped

    seq_chars = []
    stats_rows = []

    for _, row in node_df.iterrows():
        probs = {col_to_aa[c]: float(row[c]) for c in col_to_aa if c in row.index}
        if not probs:
            seq_chars.append("X")
            stats_rows.append({"site": int(row["Site"]), "map_aa": "X", "map_posterior": 0.0, "uncertain": True})
            continue
        map_aa = max(probs, key=probs.get)
        map_posterior = probs[map_aa]
        uncertain = map_posterior < posterior_threshold
        seq_chars.append(map_aa if not uncertain else "X")
        stats_rows.append(
            {
                "site": int(row["Site"]),
                "map_aa": map_aa,
                "map_posterior": map_posterior,
                "uncertain": uncertain,
            }
        )

    sequence = "".join(seq_chars)
    stats_df = pd.DataFrame(stats_rows)

    n_uncertain = stats_df["uncertain"].sum()
    n_total = len(stats_df)
    logger.info(
        "Ancestral sequence built: length %d, uncertain positions: %d (%.1f%%)",
        n_total, n_uncertain, 100.0 * n_uncertain / max(n_total, 1),
    )

    return sequence, stats_df


@click.command()
@click.option("--state-file", type=click.Path(exists=True), required=True,
              help="IQ-TREE .state file from -ancestral flag.")
@click.option("--tree-file", type=click.Path(exists=True), required=True,
              help="Rooted Newick tree (output of root_tree.py).")
@click.option("--family", required=True,
              help="Enzyme family identifier (used in FASTA header).")
@click.option("--output", "-o", type=click.Path(), required=True,
              help="Output FASTA file for ancestral sequence.")
@click.option("--output-stats", type=click.Path(), required=True,
              help="Output TSV with per-site posterior probabilities.")
@click.option("--posterior-threshold", type=float, default=0.8, show_default=True,
              help="Min posterior probability to assign MAP amino acid; below → 'X'.")
@click.option("--taxa-keywords", default=",".join(_CARNIVOROUS_TAXA_KEYWORDS),
              show_default=True,
              help="Comma-separated taxa name keywords identifying carnivorous clade.")
@click.option("--verbose", is_flag=True, help="Enable debug logging.")
def main(
    state_file: str,
    tree_file: str,
    family: str,
    output: str,
    output_stats: str,
    posterior_threshold: float,
    taxa_keywords: str,
    verbose: bool,
) -> None:
    """Extract carnivore MRCA ancestral sequence from IQ-TREE .state file."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    state_path = Path(state_file)
    tree_path = Path(tree_file)
    output_path = Path(output)
    stats_path = Path(output_stats)
    keywords = [kw.strip() for kw in taxa_keywords.split(",") if kw.strip()]

    state_df = _load_state_file(state_path)
    node_name = _find_mrca_node(tree_path, keywords)
    ancestral_seq, stats_df = _build_ancestral_sequence(state_df, node_name, posterior_threshold)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    record = SeqRecord(
        Seq(ancestral_seq),
        id=f"{family}_MRCA_ancestor",
        description=f"Ancestral MAP sequence for MRCA of carnivorous taxa | node={node_name} | threshold={posterior_threshold}",
    )
    with open(output_path, "w") as fh:
        SeqIO.write([record], fh, "fasta")
    logger.info("Wrote ancestral FASTA to %s", output_path)

    stats_path.parent.mkdir(parents=True, exist_ok=True)
    stats_df.to_csv(stats_path, sep="\t", index=False)
    logger.info("Wrote per-site stats to %s", stats_path)


if __name__ == "__main__":
    main()
