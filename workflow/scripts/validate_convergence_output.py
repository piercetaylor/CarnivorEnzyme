#!/usr/bin/env python3
"""Validate a convergent-sites TSV against the alignment and ASR that produced it.

TODO T0.1's pass criteria are deliberately not "exit code 0": the node-matching bug
class in audit/11_ancestral_reconstruction_architecture_fix.md fails *silently*,
emitting empty or fabricated rows while returning success. This script re-derives
each reported substitution from the primary files and reports what it can and
cannot confirm.

Checks:
  1. The TSV exists and has data rows beyond the header.
  2. Sampled rows trace back to real alignment columns: the reported derived_aa
     really occurs at that column in a species belonging to each reported origin,
     and the reported ancestral_aa really is a MAP state at that site in .state.
  3. Every internal node label in the ASR treefile resolves in the .state file.
  4. .asr.treefile and .asr.state agree on their node set (same IQ-TREE run).

Usage:
  python3 workflow/scripts/validate_convergence_output.py <family> [--rows N]
"""

import argparse
import random
import sys
from collections import defaultdict
from pathlib import Path

import yaml


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].strip(), []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def read_state(path):
    """{node: {site_1based: (map_aa, posterior)}} from an IQ-TREE .state file."""
    states = defaultdict(dict)
    header = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = [h.strip() for h in parts]
                i_state = header.index("State")
                i_site = header.index("Site")
                continue
            if len(parts) <= i_state:
                continue
            states[parts[0]][int(parts[i_site])] = parts[i_state]
    return states


def internal_labels(treefile):
    """Internal-node labels in a Newick file: the name following a ')'."""
    text = Path(treefile).read_text()
    labels, i = [], 0
    while True:
        i = text.find(")", i)
        if i == -1:
            break
        j = i + 1
        while j < len(text) and (text[j].isalnum() or text[j] in "._-/"):
            j += 1
        lab = text[i + 1:j]
        if lab and not lab.replace(".", "").isdigit():
            labels.append(lab)
        i = j if j > i + 1 else i + 1
    return labels


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("family")
    ap.add_argument("--rows", type=int, default=3)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    fam = args.family

    tsv = Path(f"results/convergence/{fam}.convergent_sites.tsv")
    aln = Path(f"results/alignments/{fam}.trimmed.fa")
    state = Path(f"results/phylogenies/{fam}.asr.state")
    tree = Path(f"results/phylogenies/{fam}.asr.treefile")

    failures, notes = [], []

    # --- criterion 1 -------------------------------------------------------
    print("=" * 72)
    print(f"Validating {fam}")
    print("=" * 72)
    for p in (tsv, aln, state, tree):
        print(f"  {'OK ' if p.exists() else 'MISSING'}  {p}")
        if not p.exists():
            failures.append(f"missing file: {p}")
    if failures:
        print("\nFAIL — required files absent.")
        return 1

    lines = tsv.read_text().splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, ln.split("\t"))) for ln in lines[1:] if ln.strip()]
    print(f"\n[1] TSV data rows beyond header: {len(rows)}")
    if not rows:
        failures.append("TSV has a header but zero data rows")

    alignment = read_fasta(aln)
    aln_len = len(next(iter(alignment.values()))) if alignment else 0
    print(f"    alignment: {len(alignment)} sequences x {aln_len} columns")

    states = read_state(state)
    print(f"    .state:    {len(states)} nodes")

    # --- criteria 3 & 4 ----------------------------------------------------
    labels = internal_labels(tree)
    unresolved = sorted(set(labels) - set(states))
    print(f"\n[3] internal node labels in .asr.treefile: {len(labels)}")
    print(f"    resolved against .asr.state:           {len(labels) - len(unresolved)}"
          f" / {len(labels)}")
    if unresolved:
        failures.append(f"{len(unresolved)} tree node(s) absent from .state: "
                        f"{unresolved[:5]}")
    state_only = sorted(set(states) - set(labels))
    print(f"[4] nodes in .state but not in the tree:    {len(state_only)}")
    if state_only:
        notes.append(f"nodes present only in .state: {state_only[:5]} — expected 0 if "
                     "both came from one IQ-TREE run")

    # --- criterion 2 -------------------------------------------------------
    species = yaml.safe_load(Path("config/species.yaml").read_text())
    origin_of = {}
    for entry in species.get("carnivorous_species", []) or []:
        name = entry.get("name") or entry.get("species")
        if name and entry.get("carnivory_origin") is not None:
            origin_of[name] = str(entry["carnivory_origin"])

    def leaves_for_origin(code):
        out = []
        for leaf in alignment:
            sp = leaf.split("|")[0]
            if origin_of.get(sp) == str(code):
                out.append(leaf)
        return out

    if rows:
        random.seed(args.seed)
        sample = rows if len(rows) <= args.rows else random.sample(rows, args.rows)
        print(f"\n[2] hand-checking {len(sample)} row(s) against the alignment\n")
        for r in sample:
            pos = int(r["aln_position"])
            anc, der = r["ancestral_aa"], r["derived_aa"]
            codes = [c for c in r["lineages"].split(",") if c]
            print(f"  --- position {pos}  {anc} -> {der}  "
                  f"lineages={r['lineages']}  n={r['n_lineages']}  "
                  f"cat={r.get('category','')}")
            if not (1 <= pos <= aln_len):
                failures.append(f"position {pos} outside alignment (1..{aln_len})")
                print(f"      POSITION OUT OF RANGE (alignment is {aln_len} columns)")
                continue
            col = {lf: s[pos - 1] for lf, s in alignment.items()}
            for code in codes:
                lv = leaves_for_origin(code)
                hits = [l for l in lv if col.get(l) == der]
                status = "OK " if hits else "FAIL"
                if not hits:
                    failures.append(
                        f"position {pos}: no leaf of origin {code} carries "
                        f"derived '{der}'")
                print(f"      {status} origin {code}: {len(lv)} leaf/leaves, "
                      f"{len(hits)} carrying '{der}'"
                      + (f"  e.g. {hits[0]}" if hits else ""))
                for l in lv[:4]:
                    print(f"             {col.get(l,'?')}  {l}")
            anc_nodes = [n for n, sites in states.items() if sites.get(pos) == anc]
            status = "OK " if anc_nodes else "FAIL"
            if not anc_nodes:
                failures.append(
                    f"position {pos}: ancestral '{anc}' is not a MAP state at any node")
            print(f"      {status} ancestral '{anc}' is the MAP state at "
                  f"{len(anc_nodes)} node(s)"
                  + (f", e.g. {anc_nodes[0]}" if anc_nodes else ""))
            print()

    # --- verdict -----------------------------------------------------------
    print("=" * 72)
    for n in notes:
        print(f"NOTE: {n}")
    if failures:
        print(f"FAIL — {len(failures)} problem(s):")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("PASS — sampled rows trace back to real alignment columns, and every")
    print("       internal node resolves against the .state file.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
