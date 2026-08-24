#!/usr/bin/env python3
"""Identify convergent amino acid substitutions across independent carnivorous plant lineages.

Reads a trimmed multiple sequence alignment, an IQ-TREE rooted gene tree, an IQ-TREE
.state file (ancestral sequence reconstruction), and species metadata. For each alignment
position, detects cases where ≥2 independent carnivorous lineages independently acquired
the same derived amino acid from the same ancestral state. Tests statistical significance
with a Poisson model and applies BH FDR correction.

Output TSV columns: family, aln_position, ancestral_aa, derived_aa, lineages,
n_lineages, min_posterior, p_value, q_value_bh, category
"""

import logging
import sys
from pathlib import Path
from collections import defaultdict
from itertools import combinations

import click
import numpy as np
import yaml
from scipy import stats
from Bio import SeqIO, AlignIO
import ete3

logger = logging.getLogger(__name__)

AA_COLS = list("ARNDCQEGHILKMFPSTWYV")
_AA_SET = set(AA_COLS)


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def _load_species_meta(species_yaml: Path) -> dict:
    """Return {species_code: {carnivorous, carnivory_origin, name}}."""
    with open(species_yaml, encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    meta = {}

    # Load carnivorous species
    for species_name, entry in data.get("carnivorous_species", {}).items():
        code = entry.get("code") or species_name.replace(" ", "_")
        meta[code] = {
            "carnivorous": True,
            "carnivory_origin": entry.get("carnivory_origin", 0),
            "name": species_name,
        }

    # Load outgroup (non-carnivorous) species
    for species_name, entry in data.get("outgroup_species", {}).items():
        code = entry.get("code") or species_name.replace(" ", "_")
        meta[code] = {
            "carnivorous": False,
            "carnivory_origin": 0,
            "name": species_name,
        }

    return meta


def _seq_id_to_species(seq_id: str, known_species: set) -> str | None:
    """Map a sequence ID to a species code.

    Handles 'Species_name|accession' format (pipe delimiter) and plain
    accession-only IDs. Returns None if species cannot be identified.
    """
    if "|" in seq_id:
        candidate = seq_id.split("|")[0]
        if candidate in known_species:
            return candidate
        # Try prefix matching for species with underscores
        for sp in known_species:
            if seq_id.startswith(sp + "|"):
                return sp
    # Fall back to prefix match without delimiter
    for sp in known_species:
        if seq_id.startswith(sp):
            return sp
    return None


def _parse_alignment(fasta_path: Path) -> dict[str, str]:
    """Return {seq_id: aligned_sequence_string}."""
    seqs = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        seqs[record.id] = str(record.seq)
    if not seqs:
        raise ValueError(f"No sequences found in alignment: {fasta_path}")
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) > 1:
        raise ValueError(f"Alignment has ragged lengths: {lengths}")
    return seqs


def _parse_state_file(state_path: Path) -> dict[str, dict[int, tuple[str, float]]]:
    """Parse IQ-TREE .state file.

    Returns {node_label: {site_1based: (map_aa, posterior_prob)}}.
    Skips leaf rows (only internal nodes have ancestral reconstruction).
    """
    states: dict[str, dict[int, tuple[str, float]]] = {}
    with open(state_path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#") or not line.strip():
                continue
            if header is None:
                header = line.split("\t")
                continue
            parts = line.split("\t")
            if len(parts) < 22:
                continue
            node = parts[0]
            try:
                site = int(parts[1])
            except ValueError:
                continue
            map_aa = parts[-1].strip()
            # Posterior probability of MAP state
            try:
                aa_idx = AA_COLS.index(map_aa) + 2  # columns 2..21 are p_A..p_V
                posterior = float(parts[aa_idx])
            except (ValueError, IndexError):
                posterior = 0.0
            if node not in states:
                states[node] = {}
            states[node][site] = (map_aa, posterior)
    return states


def _load_tree(tree_path: Path) -> ete3.Tree:
    """Load and return an ete3 tree, naming unnamed internal nodes."""
    tree = ete3.Tree(str(tree_path), format=1)
    counter = 1
    for node in tree.traverse("postorder"):
        if not node.is_leaf() and not node.name:
            node.name = f"Node{counter}"
            counter += 1
    return tree


# ---------------------------------------------------------------------------
# Branch length and substitution rate helpers
# ---------------------------------------------------------------------------

def _get_leaf_parent_name(tree: ete3.Tree, leaf_name: str) -> str | None:
    """Return the name of the direct parent of a leaf node."""
    results = tree.search_nodes(name=leaf_name)
    if not results:
        return None
    leaf_node = results[0]
    parent = leaf_node.up
    if parent is None:
        return None
    return parent.name


def _get_carnivorous_branch_lengths(
    tree: ete3.Tree,
    species_meta: dict,
    leaf_species_map: dict[str, str],
) -> float:
    """Sum of all branch lengths leading to carnivorous leaf nodes."""
    total = 0.0
    for leaf in tree.get_leaves():
        species = leaf_species_map.get(leaf.name)
        if species and species_meta.get(species, {}).get("carnivorous"):
            total += leaf.dist
    return total if total > 0 else 1.0


def _background_substitution_rate(alignment: dict[str, str]) -> float:
    """Estimate mean pairwise amino acid diversity across the alignment.

    Returns a per-site substitution rate (fraction of sites differing,
    averaged over all pairwise comparisons, ignoring gaps).
    """
    seqs = [s for s in alignment.values()]
    if len(seqs) < 2:
        return 0.01
    n_sites = len(seqs[0])
    pairwise_divergences = []
    for s1, s2 in combinations(seqs, 2):
        comparable = [(a, b) for a, b in zip(s1, s2) if a not in "-X" and b not in "-X"]
        if not comparable:
            continue
        diff = sum(1 for a, b in comparable if a != b)
        pairwise_divergences.append(diff / len(comparable))
    if not pairwise_divergences:
        return 0.01
    return float(np.mean(pairwise_divergences))


# ---------------------------------------------------------------------------
# Core convergence detection
# ---------------------------------------------------------------------------

def _detect_convergent_sites(
    alignment: dict[str, str],
    tree: ete3.Tree,
    states: dict[str, dict[int, tuple[str, float]]],
    species_meta: dict,
    family: str,
    threshold: float,
    min_lineages: int,
) -> list[dict]:
    """Detect convergent amino acid substitutions across carnivorous lineages.

    For each alignment column, finds positions where ≥min_lineages independent
    carnivorous origins independently acquired the same (anc_aa → derived_aa)
    substitution. Performs a Poisson significance test.
    """
    known_species = set(species_meta.keys())

    # Map leaf name → species code
    leaf_species_map: dict[str, str] = {}
    for leaf in tree.get_leaves():
        sp = _seq_id_to_species(leaf.name, known_species)
        if sp:
            leaf_species_map[leaf.name] = sp
        else:
            logger.debug("Could not map leaf '%s' to any species", leaf.name)

    carni_branch_len = _get_carnivorous_branch_lengths(tree, species_meta, leaf_species_map)
    background_rate = _background_substitution_rate(alignment)
    # Expected convergent substitutions per site under null:
    # Poisson mu = background_rate × total_carnivorous_branch_length × 1/19
    mu_per_site = background_rate * carni_branch_len / 19.0

    n_sites = len(next(iter(alignment.values())))
    results = []

    for site_0 in range(n_sites):
        site_1 = site_0 + 1  # 1-based for .state file lookup

        # Collect (carnivory_origin, anc_aa, derived_aa, posterior, leaf_name)
        # per carnivorous leaf that has a substitution at this site
        events_by_origin: dict[int, list[tuple[str, str, float, str]]] = defaultdict(list)

        for leaf in tree.get_leaves():
            species = leaf_species_map.get(leaf.name)
            if not species:
                continue
            sm = species_meta.get(species, {})
            if not sm.get("carnivorous"):
                continue
            origin = sm.get("carnivory_origin", 0)
            if origin == 0:
                continue

            derived_aa = alignment.get(leaf.name, "")[site_0] if leaf.name in alignment else None
            if derived_aa is None or derived_aa in "-X":
                continue
            if derived_aa not in _AA_SET:
                continue

            parent_name = _get_leaf_parent_name(tree, leaf.name)
            if not parent_name:
                continue

            parent_state = states.get(parent_name, {}).get(site_1)
            if parent_state is None:
                continue
            anc_aa, posterior = parent_state

            if anc_aa == derived_aa:
                continue  # no substitution
            if anc_aa not in _AA_SET:
                continue

            events_by_origin[origin].append((anc_aa, derived_aa, posterior, leaf.name))

        if not events_by_origin:
            continue

        # Deduplicate within each origin: one event per origin per (anc, derived) pair.
        # If multiple leaves in the same origin have the same substitution, keep the
        # one with the lowest posterior (most conservative).
        deduped: dict[int, dict[tuple[str, str], tuple[float, str]]] = {}
        for origin, evlist in events_by_origin.items():
            for anc_aa, derived_aa, posterior, leaf_name in evlist:
                key = (anc_aa, derived_aa)
                if origin not in deduped:
                    deduped[origin] = {}
                if key not in deduped[origin] or posterior < deduped[origin][key][0]:
                    deduped[origin][key] = (posterior, leaf_name)

        # Find (anc_aa, derived_aa) pairs spanning ≥ min_lineages origins
        pair_origins: dict[tuple[str, str], list[tuple[int, float, str]]] = defaultdict(list)
        for origin, pair_dict in deduped.items():
            for (anc_aa, derived_aa), (posterior, leaf_name) in pair_dict.items():
                pair_origins[(anc_aa, derived_aa)].append((origin, posterior, leaf_name))

        for (anc_aa, derived_aa), origin_list in pair_origins.items():
            n_origins = len(origin_list)
            if n_origins < min_lineages:
                continue

            min_posterior = min(p for _, p, _ in origin_list)

            if min_posterior < 0.5:
                logger.debug(
                    "Site %d %s→%s: min_posterior %.2f < 0.5, skipping",
                    site_1, anc_aa, derived_aa, min_posterior,
                )
                continue

            category = "full" if min_posterior >= threshold else "partial"

            lineage_codes = [str(o) for o, _, _ in sorted(origin_list)]
            lineages_str = ",".join(lineage_codes)
            leaf_names = ";".join(ln for _, _, ln in sorted(origin_list))

            # Poisson test: P(X ≥ n_origins) under null
            p_value = float(stats.poisson.sf(n_origins - 1, mu=mu_per_site))

            results.append({
                "family": family,
                "aln_position": site_1,
                "ancestral_aa": anc_aa,
                "derived_aa": derived_aa,
                "lineages": lineages_str,
                "n_lineages": n_origins,
                "min_posterior": round(min_posterior, 4),
                "p_value": p_value,
                "q_value_bh": None,  # filled in later
                "category": category,
                "_leaf_names": leaf_names,
            })

    return results


def _apply_bh_correction(results: list[dict]) -> list[dict]:
    """Apply Benjamini-Hochberg FDR correction to p_value column in-place."""
    if not results:
        return results
    p_vals = np.array([r["p_value"] for r in results])
    try:
        q_vals = stats.false_discovery_control(p_vals, method="bh")
    except AttributeError:
        # scipy < 1.11: manual BH
        n = len(p_vals)
        order = np.argsort(p_vals)
        q_vals = np.empty(n)
        q_vals[order] = p_vals[order] * n / (np.arange(n) + 1)
        # Enforce monotonicity (cumulative minimum from the right)
        q_vals = np.minimum.accumulate(q_vals[::-1])[::-1]
    for r, q in zip(results, q_vals):
        r["q_value_bh"] = round(float(q), 6)
    return results


def _write_tsv(results: list[dict], output_path: Path) -> None:
    """Write convergent sites to a tab-separated file."""
    columns = [
        "family", "aln_position", "ancestral_aa", "derived_aa",
        "lineages", "n_lineages", "min_posterior", "p_value", "q_value_bh", "category",
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in results:
            fh.write("\t".join(str(row[c]) for c in columns) + "\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@click.command()
@click.option(
    "--alignment", "-a",
    type=click.Path(exists=True),
    required=True,
    help="Trimmed multiple sequence alignment (FASTA).",
)
@click.option(
    "--tree", "-t",
    type=click.Path(exists=True),
    required=True,
    help="Rooted IQ-TREE gene tree (Newick, format 1 with branch lengths and labels).",
)
@click.option(
    "--state", "-s",
    type=click.Path(exists=True),
    required=True,
    help="IQ-TREE .state file from ancestral sequence reconstruction (--ancestral flag).",
)
@click.option(
    "--species-yaml",
    type=click.Path(exists=True),
    required=True,
    help="Path to config/species.yaml with carnivory_origin metadata.",
)
@click.option(
    "--output", "-o",
    type=click.Path(),
    required=True,
    help="Output TSV file path for convergent sites.",
)
@click.option(
    "--family",
    type=str,
    required=True,
    help="Enzyme family identifier (e.g. chitinases_gh19). Written into output TSV.",
)
@click.option(
    "--threshold",
    type=float,
    default=0.8,
    show_default=True,
    help="Minimum posterior probability to classify a site as 'full' (vs 'partial').",
)
@click.option(
    "--min-lineages",
    type=int,
    default=2,
    show_default=True,
    help="Minimum number of independent carnivorous origins required for a convergent site.",
)
@click.option(
    "--verbose",
    is_flag=True,
    help="Enable debug logging.",
)
def main(
    alignment: str,
    tree: str,
    state: str,
    species_yaml: str,
    output: str,
    family: str,
    threshold: float,
    min_lineages: int,
    verbose: bool,
) -> None:
    """Detect convergent amino acid substitutions in digestive enzyme families."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    alignment_path = Path(alignment)
    tree_path = Path(tree)
    state_path = Path(state)
    species_path = Path(species_yaml)
    output_path = Path(output)

    logger.info("Loading species metadata from %s", species_path)
    species_meta = _load_species_meta(species_path)
    n_carni = sum(1 for v in species_meta.values() if v["carnivorous"])
    logger.info("  %d species total, %d carnivorous", len(species_meta), n_carni)

    logger.info("Parsing alignment from %s", alignment_path)
    aln = _parse_alignment(alignment_path)
    n_seqs = len(aln)
    n_cols = len(next(iter(aln.values())))
    logger.info("  %d sequences × %d alignment columns", n_seqs, n_cols)

    logger.info("Loading tree from %s", tree_path)
    tree_obj = _load_tree(tree_path)
    n_leaves = len(tree_obj.get_leaves())
    logger.info("  %d leaves", n_leaves)

    logger.info("Parsing .state file from %s", state_path)
    state_data = _parse_state_file(state_path)
    n_nodes = len(state_data)
    logger.info("  %d internal nodes reconstructed", n_nodes)

    logger.info("Running convergence detection (threshold=%.2f, min_lineages=%d)",
                threshold, min_lineages)
    results = _detect_convergent_sites(
        alignment=aln,
        tree=tree_obj,
        states=state_data,
        species_meta=species_meta,
        family=family,
        threshold=threshold,
        min_lineages=min_lineages,
    )
    logger.info("  %d candidate convergent sites before FDR correction", len(results))

    results = _apply_bh_correction(results)
    n_sig = sum(1 for r in results if r["q_value_bh"] is not None and r["q_value_bh"] < 0.05)
    logger.info("  %d sites with q < 0.05 after BH correction", n_sig)

    _write_tsv(results, output_path)
    logger.info("Wrote %d convergent sites to %s", len(results), output_path)


if __name__ == "__main__":
    main()
