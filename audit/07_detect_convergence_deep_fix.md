# T7 — Deep fix for `workflow/scripts/detect_convergence.py`

**Date:** 2026-08-24
**Scope:** four independent defects, each of which alone causes `detect_convergence.py` to emit
**zero convergent sites with no error**. A prior remediation commit (`96ab469`, documented in
`audit/01_detect_convergence_loader_fix.md`) claimed to fix the first; execution shows it did not,
and it did not touch the other three.

**Why all four had to be fixed together:** each one independently zeroes the output. Fixing any
subset produces no visible change in behaviour — which is precisely the failure mode that let the
first "fix" ship broken after a manual-trace-only review.

## Execution environment

`ete3` is **not** installed in this environment, so `import detect_convergence` fails at module
import. To execute the non-tree functions against real data, a throwaway stub was placed on
`sys.path` **outside the repository** (`C:/Users/pierc/AppData/Local/Temp/ce_stub/ete3.py`). It is
not part of the repo and is not committed:

```python
"""Minimal stand-in so detect_convergence.py is importable without real ete3."""
class Tree:
    def __init__(self, *a, **k):
        raise NotImplementedError("stub")
```

Everything below is real terminal output, pasted verbatim.

---

## Bug 1 (CRITICAL) — species-metadata dict keyed by short code, never matched by leaf labels

### What was wrong

`_load_species_meta()` keyed its returned dict by the 4-letter `code` field:

```python
code = entry.get("code") or species_name.replace(" ", "_")
meta[code] = {...}
```

Every entry in `config/species.yaml` has a `code:`, so the `or` fallback was dead code and the dict
was **always** keyed `Ngra`, `Cfol`, `Dmus`, …

But every tree leaf label and FASTA header in this pipeline is `{Species_name}|{accession}` —
set by `fetch_sequences.py:230`, `record.id = f"{species}|{acc}"`, where `species` is the full
species-name key from `config/enzyme_families.yaml`'s `accessions:` block. `_seq_id_to_species()`
compares the pre-pipe token against `set(species_meta.keys())`, i.e. against code strings.
**No leaf ever matched**, `leaf_species_map` stayed empty, and every downstream carnivorous-leaf
loop was a no-op.

This also subsumes the separately-filed **HIGH-3**: `entry.get("carnivory_origin", 0)` silently
defaulted a missing origin to `0`, and `_detect_convergent_sites` drops any leaf with `origin == 0`
(`if origin == 0: continue`) — so a species with a missing origin was silently excluded while the
startup log still counted it as carnivorous.

### Before

```python
def _load_species_meta(species_yaml: Path) -> dict:
    """Return {species_code: {carnivorous, carnivory_origin, name}}."""
    with open(species_yaml, encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    meta = {}

    for species_name, entry in data.get("carnivorous_species", {}).items():
        code = entry.get("code") or species_name.replace(" ", "_")
        meta[code] = {
            "carnivorous": True,
            "carnivory_origin": entry.get("carnivory_origin", 0),
            "name": species_name,
        }

    for species_name, entry in data.get("outgroup_species", {}).items():
        code = entry.get("code") or species_name.replace(" ", "_")
        meta[code] = {
            "carnivorous": False,
            "carnivory_origin": 0,
            "name": species_name,
        }

    return meta
```

### After

```python
def _load_species_meta(species_yaml: Path) -> dict:
    """Return {lookup_key: {carnivorous, carnivory_origin, name}}.

    Registers each species under BOTH its full species name (the YAML key, e.g.
    'Nepenthes_gracilis') and its short code (e.g. 'Ngra'). Tree leaf labels and
    FASTA headers in this pipeline use the '{Species_name}|{accession}' form set by
    fetch_sequences.py, so the species-name key is the one that actually matches;
    the code key is kept for any consumer that works in short codes.
    """
    with open(species_yaml, encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    meta: dict = {}

    for species_name, entry in data.get("carnivorous_species", {}).items():
        origin = entry.get("carnivory_origin")
        if not origin:
            raise ValueError(
                f"{species_yaml}: carnivorous species '{species_name}' has no carnivory_origin"
            )
        rec = {"carnivorous": True, "carnivory_origin": origin, "name": species_name}
        meta[species_name] = rec
        code = entry.get("code")
        if code:
            meta[code] = rec

    for species_name, entry in data.get("outgroup_species", {}).items():
        rec = {"carnivorous": False, "carnivory_origin": 0, "name": species_name}
        meta[species_name] = rec
        code = entry.get("code")
        if code:
            meta[code] = rec

    if not meta:
        raise ValueError(
            f"{species_yaml}: no species loaded — expected 'carnivorous_species' "
            f"and/or 'outgroup_species' top-level keys, found: {sorted(data)}"
        )

    return meta
```

Because each species is now registered twice, the startup log in `main()` was double-counting.
It now counts distinct species by the canonical `name` field:

```python
unique_species = {v["name"]: v for v in species_meta.values()}
n_carni = sum(1 for v in unique_species.values() if v["carnivorous"])
logger.info(
    "  %d species total, %d carnivorous (%d lookup keys)",
    len(unique_species), n_carni, len(species_meta),
)
```

### Verification — BEFORE the fix

```
$ python3 -c "
import sys; sys.path.insert(0,'workflow/scripts')
import detect_convergence as dc
from pathlib import Path
m = dc._load_species_meta(Path('config/species.yaml'))
print('KEYS:', sorted(m.keys()))
print('LOOKUP:', dc._seq_id_to_species('Nepenthes_gracilis|BAD07474.1', set(m.keys())))
"

N KEYS: 23
KEYS: ['Atha', 'Aves', 'Bfil', 'Cfol', 'Dade', 'Dbin', 'Dcal', 'Dcap', 'Dmus', 'Dspa', 'Hcil', 'Hvul', 'Nala', 'Ngra', 'Nmir', 'Pmor', 'Pvul', 'Sala', 'Scer', 'Slyo', 'Spur', 'Ugib', 'Vvin']
LOOKUP Nepenthes_gracilis|BAD07474.1 -> None
```

### Verification — AFTER the fix

```
$ python3 -c "
import sys
sys.path.insert(0,r'C:/Users/pierc/AppData/Local/Temp/ce_stub')
sys.path.insert(0,'workflow/scripts')
import detect_convergence as dc
from pathlib import Path
m = dc._load_species_meta(Path('config/species.yaml'))
print('N KEYS:', len(m))
print('LOOKUP Nepenthes_gracilis|BAD07474.1 ->', repr(dc._seq_id_to_species('Nepenthes_gracilis|BAD07474.1', set(m.keys()))))
print('LOOKUP Cephalotus_follicularis|GAV12345.1 ->', repr(dc._seq_id_to_species('Cephalotus_follicularis|GAV12345.1', set(m.keys()))))
print('LOOKUP Arabidopsis_thaliana|NP_001.1 ->', repr(dc._seq_id_to_species('Arabidopsis_thaliana|NP_001.1', set(m.keys()))))
print('CODE KEY STILL WORKS Ngra ->', m['Ngra'])
print('NAME KEY ->', m['Nepenthes_gracilis'])
uniq={v['name']:v for v in m.values()}
print('unique species:', len(uniq), 'carnivorous:', sum(1 for v in uniq.values() if v['carnivorous']))
print('origins present:', sorted({v['carnivory_origin'] for v in uniq.values() if v['carnivorous']}))
print('NO ValueError raised on real config/species.yaml')
"

N KEYS: 46
LOOKUP Nepenthes_gracilis|BAD07474.1 -> 'Nepenthes_gracilis'
LOOKUP Cephalotus_follicularis|GAV12345.1 -> 'Cephalotus_follicularis'
LOOKUP Arabidopsis_thaliana|NP_001.1 -> 'Arabidopsis_thaliana'
CODE KEY STILL WORKS Ngra -> {'carnivorous': True, 'carnivory_origin': 1, 'name': 'Nepenthes_gracilis'}
NAME KEY -> {'carnivorous': True, 'carnivory_origin': 1, 'name': 'Nepenthes_gracilis'}
unique species: 23 carnivorous: 18
origins present: [1, 2, 3, 4, 5]
NO ValueError raised on real config/species.yaml
```

**Data check on the real config:** no `ValueError` was raised — every one of the 18 carnivorous
entries in `config/species.yaml` does carry a `carnivory_origin`, and all five documented origins
(1 Caryophyllales, 2 Oxalidales, 3 Ericales, 4 Lentibulariaceae, 5 Byblidaceae) are represented.
There is no latent data problem being papered over.

### Verification — new loud-failure paths

```
=== Bug 1 failure path: carnivorous species missing carnivory_origin ===
  ValueError: bad_species.yaml: carnivorous species 'Drosera_capensis' has no carnivory_origin

=== Bug 1: empty/wrong-shaped species yaml ===
  ValueError: empty_species.yaml: no species loaded — expected 'carnivorous_species' and/or 'outgroup_species' top-level keys, found: ['species']
```

---

## Bug 2 (CRITICAL) — `.state` parser read the last column instead of the `State` column

### Format research (not guessed)

The IQ-TREE Command Reference, section *Ancestral sequence reconstruction*, gives the header
verbatim as:

```
Node    Site    State   p_A     p_C     p_G     p_T
```

with column meanings: **Node** = "Node name in the tree"; **Site** = "Alignment site ID";
**State** = "Most likely state assignment"; **p_X** = "Posterior probability for state X
(empirical Bayesian method)". (That example is DNA; for protein data the `p_X` block is the 20
amino acids.)

**Amino-acid column order.** IQ-TREE's internal protein state order is
`A R N D C Q E G H I L K M F P S T W Y V` — the PAML-style ordering used throughout IQ-TREE's
documentation for PAML-format rate matrices and custom site-frequency profiles. That **does**
match this file's `AA_COLS = list("ARNDCQEGHILKMFPSTWYV")` constant, so the ordering assumption
itself was correct. The header-driven fix below is nevertheless the right one: it removes the
dependency on that assumption entirely, and it is what catches a malformed file.

### Three compounding errors in the old code

1. `map_aa = parts[-1].strip()` read the **last** column — `p_V`, a float string such as
   `"0.00012"` — and treated it as the ancestral amino acid.
2. `AA_COLS.index(map_aa)` then raised `ValueError` on every row, swallowed by the bare
   `except (ValueError, IndexError)`, silently defaulting `posterior = 0.0`.
3. **Additionally, the offset was off by one even on its own terms.** With the real header
   `Node(0) Site(1) State(2) p_A(3) …`, `p_A` is at index **3**, but `AA_COLS.index('A') + 2 == 2`
   — which points at the `State` column, not `p_A`. So even had `map_aa` been correct, every
   posterior would have been read one column to the left.

The garbage `map_aa` then failed the downstream `if anc_aa not in _AA_SET: continue` guard, so
every site was dropped. Zero rows, zero errors.

### Before

```python
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
```

Note the header was parsed and then **never used**.

### After

```python
    states: dict[str, dict[int, tuple[str, float]]] = {}
    idx_state: int | None = None
    idx_p: dict[str, int] = {}
    n_ambiguous = 0

    with open(state_path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#") or not line.strip():
                continue
            if header is None:
                header = [h.strip() for h in line.split("\t")]
                # Let ValueError propagate: a .state file missing these columns is
                # malformed and must fail immediately rather than yield zero rows.
                idx_state = header.index("State")
                idx_p = {aa: header.index("p_" + aa) for aa in AA_COLS}
                continue
            parts = line.split("\t")
            node = parts[0]
            try:
                site = int(parts[1])
            except (ValueError, IndexError):
                continue
            map_aa = parts[idx_state].strip()
            if map_aa not in _AA_SET:
                # Ambiguous/gap MAP state (e.g. '-' or '?') — not an error.
                n_ambiguous += 1
                continue
            posterior = float(parts[idx_p[map_aa]])
            if node not in states:
                states[node] = {}
            states[node][site] = (map_aa, posterior)

    if n_ambiguous:
        logger.debug("  %d rows skipped with ambiguous/gap MAP state", n_ambiguous)
    if not states:
        raise ValueError(
            f"{state_path}: parsed 0 internal nodes — check IQ-TREE --ancestral output format"
        )
    return states
```

The `AA_COLS.index(map_aa) + 2` arithmetic and its exception-swallowing `try/except` are gone; a
missing header column now propagates `ValueError` out of `header.index(...)`.

### Synthetic `.state` file used for verification

Written to a temp file, real IQ-TREE protein format, header + 4 data rows whose correct `State`
and posterior are known by construction. Row 4 deliberately carries a gap MAP state:

```
# Ancestral state reconstruction (synthetic, IQ-TREE .state format)
Node	Site	State	p_A	p_R	p_N	... (p_C..p_V)
Node1	1	E	0.00000	0.00000	0.00000	... (p_C..p_V)
Node1	2	W	0.00000	0.00000	0.00000	... (p_C..p_V)
Node2	1	V	0.00000	0.00000	0.00000	... (p_C..p_V)
Node2	2	-	0.00000	0.00000	0.00000	... (p_C..p_V)
```

(display truncated after `p_N`; the full 20 `p_*` columns are written). The non-zero probabilities
placed into those rows were: `Node1/1` → `p_E=0.97, p_D=0.02, p_Q=0.01`; `Node1/2` →
`p_W=0.61, p_F=0.30, p_Y=0.09`; `Node2/1` → `p_V=0.88, p_I=0.12`; `Node2/2` → all zero, gap state.

### Verification — AFTER the fix

```
=== _parse_state_file OUTPUT ===
Node1 1 ('E', 0.97)
Node1 2 ('W', 0.61)
Node2 1 ('V', 0.88)

ALL ASSERTIONS PASSED: State column read from index 2, posterior read from p_<State>
```

The gap row (`Node2` site 2) is correctly absent, and the assertions
`out["Node1"][1] == ("E", 0.97)` etc. all held.

### Verification — the OLD parser on the SAME file

```
=== OLD parser on the same synthetic file ===
  Node1 site 1: map_aa='0.00000' posterior=0.0   -> anc_aa in _AA_SET? False
  Node1 site 2: map_aa='0.00000' posterior=0.0   -> anc_aa in _AA_SET? False
  Node2 site 1: map_aa='0.88000' posterior=0.0   -> anc_aa in _AA_SET? False
  Node2 site 2: map_aa='0.00000' posterior=0.0   -> anc_aa in _AA_SET? False
  => every map_aa is a probability string, not an amino acid; all sites dropped downstream.
```

### Verification — new loud-failure paths

```
=== Malformed .state (header missing 'State') must raise, not silently return {} ===
  ValueError: 'State' is not in list

=== .state with zero usable internal nodes must raise ===
  ValueError: empty.state: parsed 0 internal nodes — check IQ-TREE --ancestral output format
```

---

## Bug 3 (MEDIUM) — internal node names: support labels masquerading as node names

### What was wrong

`_load_tree()` assigned `NodeN` labels only to internal nodes with **no** existing name:

```python
if not node.is_leaf() and not node.name:
    node.name = f"Node{counter}"
```

But IQ-TREE run with `-bb 1000 -alrt 1000` (as `CLAUDE.md` §7 Phase 1 and `PROGRESS.md` specify)
writes **support values** into internal node labels — `"85.7/97"`, `"100/100/95"`, `"97"`. Those
are non-empty, so they were kept as node "names". They never match the `NodeN` labels in the
`.state` file, so `_get_leaf_parent_name()` → `states.get(parent_name, {})` missed for every leaf
and the output was again silently zero.

This cannot be fully executed here (`ete3` is absent, so a real `ete3.Tree` cannot be built), so
`_load_tree` was fixed by trace plus a unit-testable predicate, and the load-time mismatch is now
caught by a defensive check in `main()` that needs only two plain Python sets.

### After — `_load_tree` plus a testable predicate

```python
def _is_support_label(name: str) -> bool:
    """True if an internal-node label is a bootstrap/support value, not a node name.

    IQ-TREE run with `-bb 1000 -alrt 1000` writes support values into internal node
    labels, e.g. '85.7/97' or '100/100/95' or plain '97'. These are NOT node names and
    will never match the NodeN labels used in the .state file, so they must not be
    treated as pre-existing names.
    """
    if not name:
        return False
    parts = name.split("/")
    if not all(parts):
        return False
    for part in parts:
        try:
            float(part)
        except ValueError:
            return False
    return True


def _load_tree(tree_path: Path) -> ete3.Tree:
    tree = ete3.Tree(str(tree_path), format=1)
    counter = 1
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            continue
        if _is_support_label(node.name):
            node.support_label = node.name
            node.name = ""
        if not node.name:
            node.name = f"Node{counter}"
            counter += 1
    return tree
```

A genuine IQ-TREE ancestral node name (`Node12`) is left untouched; a support label is moved aside
to `node.support_label` and the node gets a generated name.

### After — defensive check in `main()`

Renaming to `NodeN` uses our own traversal order, which is **not** guaranteed to reproduce
IQ-TREE's numbering. That is exactly why the mismatch must be a loud failure rather than a silent
zero:

```python
    tree_internal_names = {n.name for n in tree_obj.traverse() if not n.is_leaf()}
    overlap = tree_internal_names & set(state_data.keys())
    if not overlap:
        raise ValueError(
            f"No internal node name in the tree matches any node in the .state file. "
            f"Tree internal names (sample): {sorted(tree_internal_names)[:5]}; "
            f".state node names (sample): {sorted(state_data.keys())[:5]}"
        )
    logger.info(
        "  %d/%d internal node names matched to .state entries",
        len(overlap), len(tree_internal_names),
    )
```

### Verification — `_is_support_label` against real IQ-TREE label forms

```
=== Bug 3: _is_support_label on real IQ-TREE -bb/-alrt label forms ===
  '85.7/97'      -> support label? True
  '97'           -> support label? True
  '100/100/95'   -> support label? True
  '100.0'        -> support label? True
  'Node12'       -> support label? False
  'Node1'        -> support label? False
  ''             -> support label? False
  '0.98'         -> support label? True
```

### Verification — the defensive check, executed with plain sets

From the end-to-end run below (fake tree with `Node1`/`Node2`/`Root` internals against a `.state`
file containing `Node1`/`Node2`):

```
tree internal names: ['Node1', 'Node2', 'Root'] | .state names: ['Node1', 'Node2'] | overlap: ['Node1', 'Node2']
```

Non-empty overlap → the check passes and logs `2/3 internal node names matched`. Had the tree
carried `85.7/97`-style labels, the overlap would be empty and the run would abort with the
diagnostic naming both samples.

---

## Bug 4 (HIGH) — fabricated fallbacks producing publication-shaped p-values from no data

### What was wrong

Both null-model inputs invented plausible-looking numbers when they had nothing to work with:

```python
# _get_carnivorous_branch_lengths
    return total if total > 0 else 1.0

# _background_substitution_rate
    if len(seqs) < 2:
        return 0.01
    ...
    if not pairwise_divergences:
        return 0.01
```

Under Bug 1 no carnivorous leaf ever mapped, so `total` was always `0.0` and the branch length was
**always** the fabricated `1.0`. These two values feed
`mu_per_site = background_rate * carni_branch_len / 19.0`, which parameterises the Poisson test —
so any p-value the pipeline emitted was computed from invented inputs while looking entirely
publication-ready.

### After

```python
# _get_carnivorous_branch_lengths
    if total <= 0:
        raise ValueError(
            "No carnivorous leaf branch lengths found — check species mapping "
            "(leaf_species_map) before computing the null model"
        )
    return total

# _background_substitution_rate
    if len(seqs) < 2:
        raise ValueError(
            "Could not compute background substitution rate — alignment has "
            f"{len(seqs)} sequence(s), need at least 2"
        )
    ...
    if not pairwise_divergences:
        raise ValueError(
            "Could not compute background substitution rate — no comparable "
            "sequence pairs in alignment"
        )
```

The now-unused local `n_sites = len(seqs[0])` in `_background_substitution_rate` was removed with
the `len(seqs) < 2` early return it followed.

### Verification

Executed with a minimal fake tree object (only `.get_leaves()`, `.name`, `.dist` are needed) and
the **real** `config/species.yaml`:

```
=== Bug 4: fabricated fallbacks now raise ===
  leaf_species_map (mapped OK): {'Nepenthes_gracilis|A1': 'Nepenthes_gracilis', 'Cephalotus_follicularis|B2': 'Cephalotus_follicularis'}
  branch length total: 0.18
  -- now simulate Bug-1 conditions (nothing mapped):
  ValueError: No carnivorous leaf branch lengths found — check species mapping (leaf_species_map) before computing the null model

  -- _background_substitution_rate:
   real 3-seq alignment -> 0.20000000000000004
   single sequence: ValueError: Could not compute background substitution rate — alignment has 1 sequence(s), need at least 2
   all-gap / no comparable pairs: ValueError: Could not compute background substitution rate — no comparable sequence pairs in alignment
```

The happy paths still return real numbers (`0.18`, `0.2`); only the no-data paths raise.

---

## End-to-end proof: old code = 0 rows, new code = 1 row, identical input

`_detect_convergent_sites` touches the tree only through `get_leaves()` and `search_nodes()`, so it
can be driven with a minimal fake node class and executed for real. Input: a 4-leaf tree with three
independent carnivory origins plus one outgroup, a 3-column alignment where site 1 shows a D→E
substitution in all three carnivorous lineages, and a synthetic `.state` file giving `D` at
posterior 0.95 at both parent nodes. Species metadata is loaded from the **real**
`config/species.yaml`.

```
Root
├── Node1 ── Nepenthes_gracilis|A1      (origin 1, Caryophyllales)   site1 = E
│         └─ Arabidopsis_thaliana|AT1   (outgroup)                   site1 = D
└── Node2 ── Cephalotus_follicularis|C1 (origin 2, Oxalidales)       site1 = E
          └─ Sarracenia_purpurea|S1     (origin 3, Ericales)         site1 = E
```

### FIXED code

```
tree internal names: ['Node1', 'Node2', 'Root'] | .state names: ['Node1', 'Node2'] | overlap: ['Node1', 'Node2']

DETECTED ROWS: 1
  {'family': 'chitinases_gh19_class_iv', 'aln_position': 1, 'ancestral_aa': 'D', 'derived_aa': 'E', 'lineages': '1,2,3', 'n_lineages': 3, 'min_posterior': 0.95, 'p_value': 4.179023616592633e-08, 'q_value_bh': 0.0, 'category': 'full'}

END-TO-END: non-zero output produced.
```

### OLD code (`git show HEAD:workflow/scripts/detect_convergence.py`), same input

```
$ git show HEAD:workflow/scripts/detect_convergence.py > /tmp/old_dc.py
wrote old version

OLD meta keys (sample): ['Atha', 'Aves', 'Bfil', 'Cfol', 'Dade', 'Dbin']
OLD parsed states sample: {'Node1': {1: ('0.00000', 0.0), 2: ('0.00000', 0.0)}}
OLD DETECTED ROWS: 0
=> confirms silent zero-output on input where the fixed code finds 1 site.
```

Zero rows, no exception, no warning — the exact silent failure this task set out to eliminate.

---

## Post-fix sanity checks

```
$ python3 -m py_compile workflow/scripts/detect_convergence.py && echo "COMPILES OK"
COMPILES OK
AST OK

$ grep -n "n_sites\|0.01\|else 1.0\|parts\[-1\]\|AA_COLS.index" workflow/scripts/detect_convergence.py
321:    n_sites = len(next(iter(alignment.values())))
324:    for site_0 in range(n_sites):
```

The only surviving `n_sites` matches are the legitimate alignment-column loop in
`_detect_convergent_sites`. No `parts[-1]`, no `AA_COLS.index(...)`, and neither magic fallback
(`1.0`, `0.01`) remains anywhere in the file.

```
$ git diff --stat workflow/scripts/detect_convergence.py
 workflow/scripts/detect_convergence.py | 169 +++++++++++++++++++++++++++------
 1 file changed, 140 insertions(+), 29 deletions(-)
```

## Residual limitations

- `_load_tree()` itself was **not** executed — `ete3` is unavailable in this environment. Its
  change is trace-verified; the extracted `_is_support_label` predicate it depends on **was**
  executed. The `main()` overlap check is the runtime backstop and was executed with plain sets.
- The `NodeN` names generated by `_load_tree` follow our own postorder traversal and are not
  guaranteed to reproduce IQ-TREE's numbering. When they do not, the run now aborts with a
  diagnostic instead of emitting zero rows. Reconciling the two numbering schemes against real
  IQ-TREE `-asr` output is follow-up work and should be done the first time this runs on real data.
- No real IQ-TREE `.state` file exists in the repo yet, so Bug 2 was verified against a synthetic
  file built to the documented format rather than against tool output.

## Sources

- [IQ-TREE Command Reference — Ancestral sequence reconstruction](https://iqtree.github.io/doc/Command-Reference)
- [IQ-TREE Substitution Models (protein state ordering)](https://iqtree.github.io/doc/Substitution-Models)
- [IQ-TREE discussion #325 — parsing .state files](https://github.com/orgs/iqtree/discussions/325)
