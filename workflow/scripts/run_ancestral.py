#!/usr/bin/env python3
"""run_ancestral — Marginal ancestral state reconstruction on a FIXED, already-rooted tree.

Second IQ-TREE pass of the two-pass phylogeny architecture (see
audit/11_ancestral_reconstruction_architecture_fix.md):

  pass 1  iqtree3 -s aln -m MFP -bb 1000 -alrt 1000            -> {family}.treefile
  root    root_tree.py                                          -> {family}.rooted.plain.treefile
  pass 2  THIS SCRIPT: iqtree3 -s aln -m <model> -te <rooted>
                       -o <outgroups> --ancestral               -> {family}.asr.{treefile,state}

Why a second pass instead of `-ancestral` in pass 1:
  IQ-TREE numbers the internal nodes of the tree it *writes* (`NodeN` labels in
  `<prefix>.treefile`) and the `.state` file's `Node` column refers to THOSE labels — the
  `.state` header literally says "Ancestral state reconstruction for all nodes in
  <prefix>.treefile". Rooting the pass-1 tree afterwards produces a Newick file whose
  internal labels no longer describe the same clades, and any node introduced by the
  rooting operation has no `.state` entry at all. Running ASR *after* rooting, and then
  consuming IQ-TREE's own output treefile, removes the correspondence problem entirely:
  the tree and the `.state` file come out of the same invocation.

This script therefore also VERIFIES that correspondence before declaring success:
  * every internal node of the ASR treefile is labelled, labels are unique;
  * the label set of the ASR treefile == the node set of the .state file (exactly);
  * the unrooted split set of the ASR treefile == that of the input rooted tree
    (IQ-TREE must not have changed the topology under `-te`);
  * `-o` actually applied (IQ-TREE warns "Branch separating outgroup is not found"
    when the requested outgroup is not monophyletic in the fixed tree).

Any failure raises — no silent degradation.
"""

import logging
import re
import shutil
import subprocess
import sys
from pathlib import Path

import click
import ete3

logger = logging.getLogger(__name__)

# IQ-TREE writes exactly one of these lines into <prefix>.iqtree when ModelFinder ran.
_MODEL_LINE_RE = re.compile(
    r"^Best-fit model according to (?:BIC|AICc|AIC):\s*(\S+)\s*$"
)
# Emitted by IQ-TREE when `-o` names taxa that do not form a clade in the fixed tree.
_OUTGROUP_WARNING = "Branch separating outgroup is not found"


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def _extract_best_model(report_path: Path) -> str | None:
    """Return the model string ModelFinder selected in the pass-1 `.iqtree` report.

    Returns None when the report contains no ModelFinder verdict (i.e. pass 1 was run
    with an explicit `-m <model>` rather than `-m MFP`), in which case the caller falls
    back to the configured model.
    """
    with report_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            match = _MODEL_LINE_RE.match(line.strip())
            if match:
                return match.group(1)
    return None


def _read_outgroups(outgroups_path: Path) -> list[str]:
    """Read the comma-separated outgroup tip labels written by root_tree.py."""
    text = outgroups_path.read_text(encoding="utf-8").strip()
    if not text:
        return []
    return [t.strip() for t in text.split(",") if t.strip()]


# ---------------------------------------------------------------------------
# Tree/state correspondence verification
# ---------------------------------------------------------------------------

def _internal_labels(tree: ete3.Tree) -> list[str]:
    """Return the labels of every internal node, in preorder (root included)."""
    return [n.name for n in tree.traverse("preorder") if not n.is_leaf()]


def _unrooted_splits(tree: ete3.Tree) -> set[frozenset[str]]:
    """Return the non-trivial unrooted bipartitions of `tree` as canonical leaf sets.

    Each internal edge splits the leaves in two; the split is stored as whichever side
    is lexicographically smaller so the representation is root-independent. Trivial
    splits (a single leaf, or everything-but-one-leaf) are dropped, as is the degenerate
    split contributed by a bifurcating root. This lets a rooted tree and its unrooted
    equivalent be compared directly.
    """
    all_leaves = frozenset(leaf.name for leaf in tree.get_leaves())
    n = len(all_leaves)
    splits: set[frozenset[str]] = set()
    for node in tree.traverse("postorder"):
        if node.is_leaf() or node.is_root():
            continue
        side = frozenset(leaf.name for leaf in node.get_leaves())
        if len(side) < 2 or len(side) > n - 2:
            continue
        other = all_leaves - side
        splits.add(side if sorted(side) < sorted(other) else other)
    return splits


def _state_node_names(state_path: Path) -> set[str]:
    """Return the distinct values of the `Node` column of an IQ-TREE `.state` file."""
    names: set[str] = set()
    header_seen = False
    with state_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            if not header_seen:
                header_seen = True  # column header row
                continue
            names.add(line.split("\t", 1)[0].strip())
    return names


def _verify_correspondence(
    asr_tree_path: Path,
    state_path: Path,
    input_tree_path: Path,
) -> int:
    """Raise unless the ASR treefile and `.state` file describe exactly the same nodes.

    Returns the number of internal nodes verified.
    """
    asr_tree = ete3.Tree(str(asr_tree_path), format=1)
    input_tree = ete3.Tree(str(input_tree_path), format=1)

    labels = _internal_labels(asr_tree)
    blank = [i for i, name in enumerate(labels) if not (name or "").strip()]
    if blank:
        raise ValueError(
            f"{asr_tree_path}: {len(blank)} internal node(s) carry no label. IQ-TREE "
            f"labels every internal node when --ancestral is used, so an unlabelled node "
            f"means this is not an ancestral-reconstruction treefile."
        )
    if len(set(labels)) != len(labels):
        duplicated = sorted({name for name in labels if labels.count(name) > 1})
        raise ValueError(
            f"{asr_tree_path}: duplicated internal node labels {duplicated} — node "
            f"identity is ambiguous and .state lookups cannot be trusted."
        )

    state_nodes = _state_node_names(state_path)
    tree_nodes = set(labels)
    missing_in_state = sorted(tree_nodes - state_nodes)
    missing_in_tree = sorted(state_nodes - tree_nodes)
    if missing_in_state or missing_in_tree:
        raise ValueError(
            f"Node sets of {asr_tree_path} and {state_path} disagree. "
            f"In tree but not .state: {missing_in_state[:10]}; "
            f"in .state but not tree: {missing_in_tree[:10]}."
        )

    asr_splits = _unrooted_splits(asr_tree)
    input_splits = _unrooted_splits(input_tree)
    if asr_splits != input_splits:
        only_asr = len(asr_splits - input_splits)
        only_input = len(input_splits - asr_splits)
        raise ValueError(
            f"IQ-TREE changed the topology despite -te: {only_input} split(s) of "
            f"{input_tree_path} are absent from {asr_tree_path} and {only_asr} split(s) "
            f"are new. The fixed-tree pass must preserve the rooted tree's topology."
        )

    return len(labels)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@click.command()
@click.option("--alignment", "-a", type=click.Path(exists=True), required=True,
              help="Trimmed multiple sequence alignment (FASTA) — same file used in pass 1.")
@click.option("--tree", "-t", type=click.Path(exists=True), required=True,
              help="Rooted tree with internal labels stripped (root_tree.py --output-plain).")
@click.option("--report", type=click.Path(exists=True), required=True,
              help="Pass-1 IQ-TREE report ({family}.iqtree); the ModelFinder verdict is read from it.")
@click.option("--outgroups", type=click.Path(exists=True), default=None,
              help="File of comma-separated outgroup tip labels (root_tree.py --output-outgroups). "
                   "Passed to IQ-TREE as -o so the ASR treefile is written outgroup-first.")
@click.option("--prefix", type=str, required=True,
              help="IQ-TREE output prefix, e.g. results/phylogenies/rnase_t2.asr")
@click.option("--fallback-model", type=str, default="LG+G4", show_default=True,
              help="Model to use when the pass-1 report has no ModelFinder verdict.")
@click.option("--iqtree-binary", type=str, default="iqtree3", show_default=True,
              help="IQ-TREE executable name or path.")
@click.option("--threads", type=int, default=1, show_default=True,
              help="Threads passed to IQ-TREE (-nt).")
@click.option("--seed", type=int, default=42, show_default=True,
              help="Random seed passed to IQ-TREE (-seed).")
@click.option("--allow-outgroup-nonmonophyly", is_flag=True,
              help="Downgrade the 'outgroup is not monophyletic' failure to a warning. "
                   "The ASR treefile is then written in an arbitrary orientation and the "
                   "carnivore MRCA extracted downstream may be wrong.")
@click.option("--verbose", is_flag=True, help="Enable debug logging.")
def main(
    alignment: str,
    tree: str,
    report: str,
    outgroups: str | None,
    prefix: str,
    fallback_model: str,
    iqtree_binary: str,
    threads: int,
    seed: int,
    allow_outgroup_nonmonophyly: bool,
    verbose: bool,
) -> None:
    """Run IQ-TREE ancestral state reconstruction on a fixed rooted topology."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    alignment_path = Path(alignment)
    tree_path = Path(tree)
    report_path = Path(report)
    prefix_path = Path(prefix)
    prefix_path.parent.mkdir(parents=True, exist_ok=True)

    if shutil.which(iqtree_binary) is None and not Path(iqtree_binary).exists():
        raise FileNotFoundError(
            f"IQ-TREE binary '{iqtree_binary}' not found on PATH. Install it "
            f"(conda: iqtree>=3.0) or pass --iqtree-binary with an explicit path."
        )

    model = _extract_best_model(report_path)
    if model is None:
        model = fallback_model
        logger.warning(
            "No ModelFinder verdict in %s — falling back to model '%s'. Pass 1 was "
            "probably run with an explicit -m instead of -m MFP.",
            report_path, model,
        )
    else:
        logger.info("Reusing pass-1 ModelFinder verdict: %s", model)

    cmd = [
        iqtree_binary,
        "-s", str(alignment_path),
        "-m", model,
        "-te", str(tree_path),          # fixed topology: no tree search is performed
        "--ancestral",                  # marginal ASR by empirical Bayes
        "-nt", str(threads),
        "-seed", str(seed),
        "-redo",
        "-pre", str(prefix_path),
    ]
    # NOTE: -bb/-B (UFBoot) is deliberately absent. IQ-TREE rejects it outright here:
    # "Ultrafast bootstrap does not work with -fast, -te or -n option". Support values
    # belong to the pass-1 topology search and are carried by {family}.treefile.

    outgroup_tips: list[str] = []
    if outgroups:
        outgroup_tips = _read_outgroups(Path(outgroups))
        if outgroup_tips:
            cmd += ["-o", ",".join(outgroup_tips)]
            logger.info("Orienting output tree on %d outgroup tip(s)", len(outgroup_tips))
        else:
            logger.warning("%s is empty — running without -o", outgroups)

    logger.info("Running: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        logger.error("IQ-TREE stdout tail:\n%s", proc.stdout[-4000:])
        logger.error("IQ-TREE stderr tail:\n%s", proc.stderr[-4000:])
        raise RuntimeError(
            f"IQ-TREE exited with status {proc.returncode} (prefix {prefix_path})"
        )

    asr_tree_path = Path(f"{prefix_path}.treefile")
    state_path = Path(f"{prefix_path}.state")
    for path in (asr_tree_path, state_path):
        if not path.exists():
            raise FileNotFoundError(
                f"IQ-TREE reported success but {path} was not written — check "
                f"{prefix_path}.log"
            )

    if outgroup_tips and _OUTGROUP_WARNING in proc.stdout:
        message = (
            f"IQ-TREE reported '{_OUTGROUP_WARNING}': the outgroup tips "
            f"{outgroup_tips[:5]} do not form a clade in {tree_path}, so {asr_tree_path} "
            f"was NOT written outgroup-first. Ancestral states themselves are unaffected "
            f"(marginal ASR under a reversible model is root-independent), but any "
            f"downstream MRCA lookup on this tree is unreliable."
        )
        if allow_outgroup_nonmonophyly:
            logger.warning(message)
        else:
            raise ValueError(
                message + " Re-run with --allow-outgroup-nonmonophyly to proceed anyway."
            )

    n_nodes = _verify_correspondence(asr_tree_path, state_path, tree_path)
    logger.info(
        "Verified: %d internal nodes labelled in %s and present in %s; topology matches %s",
        n_nodes, asr_tree_path, state_path, tree_path,
    )
    logger.info("Ancestral reconstruction complete (prefix %s)", prefix_path)


if __name__ == "__main__":
    main()
