# T11 — Architectural fix: ancestral reconstruction now runs on the tree the pipeline consumes

**Date:** 2026-08-24
**Files changed:** `workflow/rules/phylogeny.smk`, `workflow/rules/convergence.smk`,
`workflow/rules/ancestral_structure.smk`, `workflow/scripts/detect_convergence.py`,
`workflow/scripts/root_tree.py`, `workflow/scripts/run_ancestral.py` (new), `Snakefile`,
`config/config.yaml`

---

## 0. Verification status up front

This project has a history of "verified" claims that failed when someone executed the code, so
this section states exactly what was and was not run.

| Claim | How it was established |
|---|---|
| `-te` + `--ancestral` is a supported IQ-TREE workflow | Official docs (quoted §1) **and** executed |
| `.state` `Node` column names the nodes of `<prefix>.treefile` | `.state` header text **and** executed |
| Feeding a rooted tree to `-te` silently unroots it (reversible model) | IQ-TREE source, lead-author forum post, **and** executed |
| Rerooting after ASR breaks the tree↔`.state` correspondence | **Executed** — reproduced the shipped bug with real IQ-TREE output |
| `-B`/`-bb` is rejected together with `-te` | **Executed** — exact error message captured |
| `-alrt` + `--ancestral` drops the root node's label from `.treefile` | **Executed** — reproduced independently, twice |
| Marginal ASR is root-placement-invariant here | **Executed** — max \|Δp\| = 1e-5 across all 7 nodes × 300 sites |
| The new two-pass chain produces a self-consistent tree + `.state` | **Executed** end-to-end with real IQ-TREE 3.1.3 |
| New guard rails raise on wrong inputs | **Executed** — 4 negative tests, all raise |
| A planted 3-origin convergent substitution is recovered | **Executed** — positive control, exact match |
| Snakemake DAG builds with correct rule chaining, no ambiguity | **Executed** — `snakemake -n` dry run |
| Snakemake *executes* the new rules end-to-end | **NOT verified.** Snakemake's shell layer is broken on this Windows laptop for even a trivial `echo` rule (see §7). Every shell command in the new rules was executed manually instead. |
| Behaviour with the project's real alignments (100+ seqs, 1000+ columns) | **NOT verified.** No real alignment exists in the repo yet. |

**IQ-TREE binary actually used:** `iqtree-3.1.3-Windows.zip`, downloaded from
<https://github.com/iqtree/iqtree3/releases/download/v3.1.3/iqtree-3.1.3-Windows.zip>
(10,672,551 bytes), extracted and run standalone.

```
$ ./iqtree-3.1.3-Windows/bin/iqtree3.exe --version
IQ-TREE version 3.1.3 for Windows 64-bit built Jun 19 2026
```

Test data: a 9-taxon, 300-column amino-acid alignment simulated with IQ-TREE's own AliSim
(`--alisim sim -t true_tree.nwk -m LG --length 300 --seqtype AA --seed 7`) along a tree whose tip
labels use this pipeline's real `{Species_name}|{accession}` convention and cover carnivory
origins 1, 2, 3 and 4 plus the three configured outgroups. Origin 2 is represented by
`Cephalotus_follicularis` alone and origin 3 by `Sarracenia_purpurea` alone, matching the
single-species-origin fragility in `config/species.yaml`.

---

## 1. Step 1 — What IQ-TREE actually documents

> Note: `iqtree.org` currently serves an **expired TLS certificate**, so `http://www.iqtree.org/doc/*`
> fetches fail. The identical official documentation is served at `https://iqtree.github.io/doc/*`;
> all URLs below use that mirror.

### 1.1 `-t` vs `-te`

<https://iqtree.github.io/doc/Command-Reference> — *General options*:

> `-t` — "Specify a file containing starting tree for tree search. The special option `-t BIONJ`
> starts tree search from BIONJ tree and `-t RANDOM` starts tree search from completely random
> tree. DEFAULT: 100 parsimony trees + BIONJ tree"

> `-te` — "Like `-t` but fixing user tree. That means, no tree search is performed and IQ-TREE
> computes the log-likelihood of the fixed user tree."

Same page, *Tree topology tests*:

> `-te` — "Specify a fixed user tree to estimate model parameters. Thus it behaves like `-n 0` but
> uses a user-defined tree instead of parsimony tree."

The IQ-TREE 3 built-in help exposes the long form:

```
  --tree-fix           Fix -t tree (no tree search performed)
```

`-te` is **not** `-g`: `-g` is documented separately as a *topological constraint* tree that "can
be a multifurcating tree and need not to include all taxa". We want `-te` (a fixed tree), not `-g`.

### 1.2 `-te` combined with ancestral reconstruction — explicitly documented

<https://iqtree.github.io/doc/Command-Reference> — *Ancestral sequence reconstruction*:

> "This feature is newly introduced in version 1.6. You can combine this feature with `-te` option
> to determine ancestral sequences along a user-defined tree (Otherwise, IQ-TREE computes ancestral
> sequences of the ML tree)."

And, critically for this fix, the `-te` entry **inside that same ASR section**:

> `-te` — "Specify a user-defined tree to determine ancestral sequences along this tree. You can
> assign each node of this tree with a node name, and IQ-TREE will report the ancestral sequences
> of the corresponding nodes. If nodes do not have names, IQ-TREE will automatically assign node
> names as Node1, Node2, etc."

The feature originates in a request written by IQ-TREE's lead author, Bui Quang Minh
(<https://github.com/Cibiv/IQ-TREE/issues/9>, 2016-09-08), whose example command is exactly the
pattern adopted here:

> ```
> iqtree -s alignment -te user_tree -asr
> ```

A published third-party tool, `ancseq` (<https://github.com/YuSugihara/ancseq>), ships this same
invocation verbatim in `ancseq/ancseq.py`:

```python
cmd = (
    f"iqtree -asr "
    f"-s {seq} "
    f"-te {tree} "
    ...
    f"-m {self.args.model} "
)
```

`-asr` and `--ancestral` are the same flag (parser: `strcmp(argv[cnt], "-asr") == 0 ||
strcmp(argv[cnt], "--ancestral") == 0`), identical in IQ-TREE 2 and 3.

### 1.3 Model specification with `-te`

`-m` is **not** required — ModelFinder runs by default (Tutorial: *"Starting with version 1.5.4,
`-m MFP` is the default behavior"*), and it will use the fixed tree. **But for a reproducible
pipeline `-m` must be passed explicitly**, otherwise pass 2 can select a different model than
pass 1 did. Our `run_ancestral.py` therefore parses pass 1's ModelFinder verdict out of
`{family}.iqtree` and passes it to `-m`.

Branch lengths **are re-optimised by default** on the fixed topology; `-blfix` turns that off:

> `-blfix` — "Fix branch lengths of tree passed via `-t` or `-te`. This is useful to evaluate the
> log-likelihood of an input tree with fixed tolopogy and branch lengths. DEFAULT: OFF"

We deliberately leave re-optimisation ON: `root_tree.py` (Bio.Phylo `root_with_outgroup`) can emit
a degenerate zero-length branch at the new root, and re-optimisation repairs it. There is no
documentation page stating that `-te` re-optimises branch lengths — it is only implied by
`-blfix`'s "DEFAULT: OFF" plus the FAQ's consensus-tree example, and was confirmed by execution.

### 1.4 `.state` node naming and the file it refers to

Verbatim `.state` header from our own run (`p3.state`), matching the documented example:

```
# Ancestral state reconstruction for all nodes in .tmp_iqtree/p3.treefile
# This file can be read in MS Excel or in R with command:
#   tab=read.table('.tmp_iqtree/p3.state',header=TRUE)
# Columns are tab-separated with following meaning:
#   Node:  Node name in the tree
#   Site:  Alignment site ID
#   State: Most likely state assignment
#   p_X:   Posterior probability for state X (empirical Bayesian method)
Node	Site	State	p_A	p_R	p_N	...	p_V
```

**This header line is the entire documented specification of the correspondence**, and it names
`<prefix>.treefile` — the tree IQ-TREE itself writes — and nothing else. There is no prose
documentation page anywhere on iqtree.github.io stating that `<prefix>.treefile` is rewritten with
`NodeN` labels. Source-level confirmation (`main/treetesting.cpp`, `printAncestralSequences()`,
identical in iqtree2 and iqtree3 `master`):

```cpp
out << "# Ancestral state reconstruction for all nodes in " << tree->params->out_prefix << ".treefile" << endl
...
    // set node name if neccessary
    if (node->name.empty() || !isalpha(node->name[0])) {
        node->name = "Node" + convertIntToString(node->id-tree->leafNum+1);
```

This runs inside `printMiscInfo()` (phyloanalysis.cpp:3128) **before** `printResultTree()`
(phyloanalysis.cpp:3235), which is why the labels end up in `.treefile`.

**Consequence, and the design pivot of this fix:** the node numbering is deterministic and derived
from IQ-TREE's own output tree, and — decisively — *we never have to reproduce that numbering*,
because IQ-TREE writes the labels into a Newick file we can read directly.

### 1.5 Rooted input to `-te`

There is no Command Reference statement about rooted `-te` input. The authoritative answer is from
Minh Bui on the IQ-TREE Google Group
(<https://groups.google.com/g/iqtree/c/jOn0BUN2dGI>, 2021-06-10):

> "Assuming so: if the model is reversible, then IQ-TREE will automatically convert it to an
> unrooted tree. As the root is removed, it's not possible to obtain the root sequence."

Source (`main/treetesting.cpp`):

```cpp
if (tree->rooted && tree->getModelFactory()->isReversible()) {
    if (tree->leafNum != tree->aln->getNSeq()+1)
        outError("Tree does not have same number of taxa as alignment");
    tree->convertToUnrooted();
}
```

and rooted-detection at `tree/mtree.cpp:742`: `// 2018-01-05: assuming rooted tree if root node has two children`.

So: **a bifurcating rooted tree handed to `-te` under LG/WAG/JTT is accepted, logged as
`rooted tree`, then silently unrooted, and the root vertex introduced by rooting is dissolved.**
Verified by execution — see §3.4. There is **no** message reading "Rooted tree has been unrooted";
that string does not exist in IQ-TREE.

This does not harm ancestral *state assignment* at any internal vertex (§3.5 shows the states are
numerically identical), but it does mean the rooted tree's node structure must never be assumed to
survive into `.state`.

Also worth recording, since Phase 4A reconstructs an ancestor: Georg Hochberg, on that same thread,
warns against reconstructing the **root** node specifically ("The position of the root node along
the root branch is unknown ... The community is becoming intolerant of this mistake"). This project
is safe on that point — the carnivore MRCA is an internal node, not the root — but if anyone later
wants the true root sequence, the documented workarounds are Minh Bui's `-keep_empty_seq`
virtual-root trick or a non-reversible model (`--model-joint NONREV`; see Naser-Khdour et al. 2021,
<https://doi.org/10.1093/sysbio/syab067>).

### 1.6 `-bb`/`-alrt` in the fixed-tree pass

`-B`/`-bb` is **hard-blocked** with `-te`. Executed on 3.1.3:

```
$ iqtree3 -s sim.fa -m LG -te rooted_clean.nwk --ancestral -B 1000 -alrt 1000 ...
Ultrafast bootstrap does not work with -fast, -te or -n option
```

Source (`utils/tools.cpp`): `throw("Ultrafast bootstrap does not work with -te option");`. This is
**not documented anywhere on the web** — source and runtime only.

`-alrt` *is* allowed with `-te` (no guard in the parser). **But it must not be used in the ASR
pass**, for a reason found by execution: `-alrt` together with `--ancestral` joins the node name to
the support values (`Node2/98.6/100`) **and drops the outermost node's label entirely** while
`.state` still contains its rows. Our own pass-1 output shows exactly this:

```
...)Node4/100/100:0.2464457461)Node2/98.6/100:0.1291567695);      <- no root label
$ awk 'NR>1 && $0 !~ /^#/ {print $1}' p1.state | sort -u
Node Node1 Node2 Node3 Node4 Node5 Node6 Node7                    <- Node1 exists in .state
```

Support values belong to the topology search anyway, so the ASR pass omits both flags.

---

## 2. Step 2 — Execution was achieved

Yes. Full real execution against IQ-TREE 3.1.3 on this machine; no doc-only inference was needed
for any load-bearing claim. `which iqtree iqtree2 iqtree3` found nothing initially, so the official
Windows release binary was downloaded and run standalone (no installer, no DLL setup beyond the
bundled `libiomp5md.dll`).

Two local environment notes:

* `ete3` 3.1.3 is **broken on Python 3.13** (`ModuleNotFoundError: No module named 'cgi'` — `cgi`
  was removed from the stdlib in 3.13). Installing `legacy-cgi` restores it. This is a local
  laptop-only issue: `environment.yml` pins `python=3.11`, where `cgi` still exists. Worth knowing
  if anyone tries to run these scripts outside the conda env.
* Snakemake's shell execution is broken on this machine (§7), so DAG validation is by dry run and
  every rule's shell body was executed by hand.

---

## 3. The bug, reproduced with real IQ-TREE output

### 3.1 Pass 1 as previously shipped

`iqtree3 -s sim.fa -m LG -B 1000 -alrt 1000 -nt 2 -seed 42 --ancestral -redo -pre p1`

```
p1.treefile:
(Nepenthes_gracilis|N1:0.24,Drosera_capensis|D1:0.19,((Cephalotus_follicularis|C1:0.34,
Sarracenia_purpurea|S1:0.32)Node3/99.7/100:0.19,((Utricularia_gibba|U1:0.42,
Pinguicula_vulgaris|P1:0.27)Node5/99.6/100:0.20,((Arabidopsis_thaliana|A1:0.33,
Vitis_vinifera|V1:0.23)Node7/99.7/100:0.21,Solanum_lycopersicum|L1:0.39)Node6/99.5/100:0.18)
Node4/100/100:0.25)Node2/98.6/100:0.13);
```

Note the label form: **`Node3/99.7/100`**, i.e. `NodeN` **joined to** the aLRT/UFBoot values —
*not* a bare support label. `_is_support_label("Node3/99.7/100")` splits on `/`, tries
`float("Node3")`, fails, and returns `False`. So the existing loader kept `Node3/99.7/100` as the
node's name, which matches nothing in `.state` (whose column reads `Node3`).

### 3.2 Rooting slides the labels onto different clades

`python3 workflow/scripts/root_tree.py --treefile p1.treefile ... --output p1.rooted.treefile`

```
((((Nepenthes_gracilis|N1,Drosera_capensis|D1),(Cephalotus_follicularis|C1,
Sarracenia_purpurea|S1)Node3/99.7/100)Node2/98.6/100,(Utricularia_gibba|U1,
Pinguicula_vulgaris|P1)Node5/99.6/100)Node4/100/100,(Arabidopsis_thaliana|A1,
Vitis_vinifera|V1)Node7/99.7/100,Solanum_lycopersicum|L1)Node6/99.5/100;
```

| label | clade in `p1.treefile` (what `.state` describes) | clade in `p1.rooted.treefile` |
|---|---|---|
| `Node2` | (C,S),((U,P),((A,V),L)) — 7 taxa | (N,D),(C,S) — 4 taxa |
| `Node4` | (U,P),((A,V),L) — 5 taxa | ((N,D),(C,S)),(U,P) — 6 taxa |
| `Node6` | (A,V),L — 3 taxa | the root |
| — | *(root, unlabelled in Newick, = `Node1` in `.state`)* | (N,D) — **now unlabelled** |

### 3.3 The overlap heuristic passes on one coincidental match

Executing the shipped code against these real files:

```
$ python3 -c "... dc._load_tree('p1.rooted.treefile') ... dc._parse_state_file('p1.state') ..."
tree internal names : ['Node1', 'Node2/98.6/100', 'Node3/99.7/100', 'Node4/100/100',
                       'Node5/99.6/100', 'Node6/99.5/100', 'Node7/99.7/100']
state node names    : ['Node1', 'Node2', 'Node3', 'Node4', 'Node5', 'Node6', 'Node7']
OVERLAP (buggy check): ['Node1']

 leaf Nepenthes_gracilis|N1        -> parent Node1
 leaf Drosera_capensis|D1          -> parent Node1
 leaf Cephalotus_follicularis|C1   -> parent Node3/99.7/100
 leaf Sarracenia_purpurea|S1       -> parent Node3/99.7/100
 leaf Utricularia_gibba|U1         -> parent Node5/99.6/100
 leaf Pinguicula_vulgaris|P1       -> parent Node5/99.6/100
 leaf Arabidopsis_thaliana|A1      -> parent Node7/99.7/100
 leaf Vitis_vinifera|V1            -> parent Node7/99.7/100
 leaf Solanum_lycopersicum|L1      -> parent Node6/99.5/100
```

The single `Node1` match comes from `_load_tree`'s own postorder counter naming the newly-unlabelled
`(N,D)` node `Node1`. The overlap check passes on that one hit while **6 of 7 nodes resolve to
nothing**, and the full CLI then exits 0, having reported it as a non-alarming INFO line:

```
INFO __main__:   7 internal nodes reconstructed
INFO __main__:   1/7 internal node names matched to .state entries
INFO __main__:   0 candidate convergent sites before FDR correction
INFO __main__: Wrote 0 convergent sites to BUGGY.tsv
```

Seven of nine leaves — including every leaf carrying carnivory origins 2, 3 and 4 — contribute
nothing, and the run looks clean.

That the surviving `Node1` happened to be the *right* vertex here is luck, not correctness: it
depends entirely on which internal node the postorder counter reaches first. Nothing in the code
enforces it, and with a bifurcating (single-outgroup) rooting the newly created root vertex has no
`.state` entry at all while still being handed a `NodeN` name that will resolve to some unrelated
real node.

### 3.4 Rooted trees are silently unrooted by `-te`

Rooting on a single outgroup tip (`Solanum_lycopersicum|L1`) gives a genuine bifurcating root
(8 internal vertices). Feeding it to `-te`:

```
$ iqtree3 -s sim.fa -m LG -te single_rooted_clean.nwk --ancestral ... -pre p6
Reading input tree file single_rooted_clean.nwk ... rooted tree

p6.treefile: (Nepenthes_gracilis|N1:0.24,Drosera_capensis|D1:0.19,(...)Node2:0.13)Node1;
p6.state nodes: Node1 Node2 Node3 Node4 Node5 Node6 Node7      <- 7, not 8
```

IQ-TREE detects the rooted tree, unroots it, dissolves the root vertex, and writes an unrooted
output. This is decisive: **the rooted tree's node structure can never be assumed to match
`.state`; only IQ-TREE's own ASR output tree can.**

### 3.5 Marginal ancestral states are root-placement-invariant

Comparing the pass-1 unrooted run (`p1.state`) against the `-te`-on-rooted-tree run (`p2.state`),
vertex by vertex, across all 300 sites:

```
  Node1  vs Node1                identical MAP at 300/300 sites, max |dP| = 1e-05
  Node2  vs Node2/98.6/100       identical MAP at 300/300 sites, max |dP| = 1e-05
  Node3  vs Node3/99.7/100       identical MAP at 300/300 sites, max |dP| = 1e-05
  Node4  vs Node4/100/100        identical MAP at 300/300 sites, max |dP| = 1e-05
  Node5  vs Node5/99.6/100       identical MAP at 300/300 sites, max |dP| = 1e-05
  Node6  vs Node6/99.5/100       identical MAP at 300/300 sites, max |dP| = 1e-05
  Node7  vs Node7/99.7/100       identical MAP at 300/300 sites, max |dP| = 1e-05
```

1e-05 is the `.state` file's own printing precision. This is the expected behaviour of marginal ASR
under a reversible model, and it matters for two reasons: the fix does not change any ancestral
state *value* (it changes which node each state is attributed to), and `detect_convergence.py`'s
leaf-vs-adjacent-vertex comparison is therefore unaffected by IQ-TREE's choice of output root.

*(This experiment also showed IQ-TREE copying the input tree's internal labels straight through into
`.state` — `Node2/98.6/100` appears as a `.state` node name. That is why `root_tree.py` now strips
internal labels before the ASR pass: otherwise `.state` inherits support values that describe the
wrong clade.)*

---

## 4. The fix

### 4.1 Architecture

```
  results/alignments/{family}.trimmed.fa
              │
   [infer_tree]  iqtree3 -m MFP -bb 1000 -alrt 1000        (NO --ancestral)
              ├─> {family}.treefile      support labels, topology search result
              └─> {family}.iqtree        ModelFinder verdict
              │
   [root_tree]   root_tree.py
              ├─> {family}.rooted.treefile        rooted, labelled (inspection/figures)
              ├─> {family}.rooted.plain.treefile  rooted, labels stripped  ──┐
              └─> {family}.outgroups.txt          matched outgroup tips   ──┤
              │                                                              │
   [ancestral_reconstruction]  run_ancestral.py  <──────────────────────────┘
              │   iqtree3 -m <pass-1 model> -te <rooted.plain> -o <outgroups> --ancestral
              ├─> {family}.asr.treefile   NodeN labels, written BY THIS RUN
              └─> {family}.asr.state      Node column == those labels
              │
   [detect_convergence]  consumes .asr.treefile + .asr.state  (one consistent pair)
   [extract_ancestor]    consumes .asr.treefile + .asr.state  (one consistent pair)
```

There is no longer any node-correspondence problem to solve: the tree and the `.state` file come
out of the same IQ-TREE invocation, and the labels are read out of the tree rather than
reconstructed.

### 4.2 Why `-o` is passed

IQ-TREE writes its ASR tree unrooted. `-o` ("Specify an outgroup taxon name to root the tree. The
output tree in `.treefile` will be rooted accordingly") re-orients that output so the carnivore
clade is a single subtree, which is what `extract_ancestor.py`'s `get_common_ancestor()` needs.
Executed: with `-o`, the ASR tree came out as
`((((N,D)Node5,(C,S)Node6)Node4,(U,P)Node7)Node3,(A,V)Node1,L)Node2;` — `Node3` is exactly the
carnivore MRCA. Note that `-o` *changes* the node numbering (compare the no-`-o` run, where the
same clade is unlabelled and the root is `Node1`), which is further evidence that numbering must
never be assumed and must always be read from the output tree.

If the outgroup is not monophyletic in the fixed tree, IQ-TREE warns
`WARNING: Branch separating outgroup is not found` and falls back to its default orientation
(executed, non-fatal, exit 0). `run_ancestral.py` detects that string and **raises**, because the
downstream MRCA lookup would be unreliable; `--allow-outgroup-nonmonophyly` downgrades it to a
warning.

### 4.3 `detect_convergence.py`

* `_is_support_label()` and the old `_load_tree()` counter are **deleted**. Nothing in this script
  invents a node name any more — the previous `NodeN` counter shared a vocabulary with IQ-TREE's
  own unrelated `NodeN` numbering, which is what made a coincidental match possible.
* New `_load_asr_tree()` **raises** if any internal node is unlabelled or any label is duplicated.
  This is what catches every wrong-tree input, including IQ-TREE's own undocumented
  `-alrt + --ancestral` missing-root-label case.
* New `_verify_tree_state_correspondence()` replaces the set-overlap heuristic with a node-by-node
  check: every internal node must resolve in `.state`, **and** every leaf's direct parent must
  resolve. Both raise.
* New `_leaf_parent_names()` computes leaf→parent once (also removing an `O(n)` `search_nodes()`
  call from the inner per-site loop) and raises on duplicate leaf labels.
* Every tree leaf must be present in the alignment — previously a missing leaf was silently skipped
  by `alignment.get(leaf.name, "")`.
* **Unmapped-leaf visibility (the additional fix requested):** `_map_leaves_to_species()` now
  returns the unmapped leaves and they are logged at `WARNING` with their names, instead of one
  `debug` line per leaf.
* **Origin-coverage check:** `_origins_represented()` logs one INFO line per carnivory origin with
  its leaf count, then the run **raises** if fewer distinct origins are represented than
  `--min-lineages`. This matters because `config/species.yaml` gives origin 2 exactly one species
  (`Cephalotus_follicularis`) and origin 5 exactly one (`Byblis_filifolia`) — a single unparsed
  leaf label silently removes a whole independent origin. `--allow-insufficient-origins` downgrades
  it to `ERROR` for the genuinely single-origin `methods_benchmark` families.

All four tier1 families carry ≥4 origins (`chitinases_gh19_class_iv`: 1,2,3,4,5;
`purple_acid_phosphatase`, `rnase_t2`, `aspartic_proteases_a1b_homology`: 1,2,3,4), verified by
parsing the real configs, so raising by default is safe for the pipeline as configured.

### 4.4 `Snakefile` — wildcard constraints

`{family}` is now pinned to the tier1 keys. Without it, `results/phylogenies/rnase_t2.asr.treefile`
is also matchable by `infer_tree` with `family="rnase_t2.asr"` (and `.rooted.treefile` with
`family="rnase_t2.rooted"` — an ambiguity that already existed before this change).

---

## 5. End-to-end execution of the fixed chain

```
$ iqtree3 -s aln.fa -m MFP -bb 1000 -alrt 1000 -nt 2 -seed 42 -redo -pre fam
$ python3 workflow/scripts/root_tree.py -t fam.treefile -s config/species.yaml \
      -f config/enzyme_families.yaml --family rnase_t2 -o fam.rooted.treefile \
      --output-plain fam.rooted.plain.treefile --output-outgroups fam.outgroups.txt
INFO: Outgroup tips found in tree (3): ['Arabidopsis_thaliana|A1', 'Vitis_vinifera|V1', 'Solanum_lycopersicum|L1']
INFO: Wrote label-stripped rooted tree to fam.rooted.plain.treefile
INFO: Wrote 3 outgroup tip label(s) to fam.outgroups.txt

$ python3 workflow/scripts/run_ancestral.py -a aln.fa -t fam.rooted.plain.treefile \
      --report fam.iqtree --outgroups fam.outgroups.txt --prefix fam.asr --threads 2 --verbose
INFO: Reusing pass-1 ModelFinder verdict: LG
INFO: Orienting output tree on 3 outgroup tip(s)
INFO: Running: iqtree3 -s aln.fa -m LG -te fam.rooted.plain.treefile --ancestral -nt 2
      -seed 42 -redo -pre fam.asr -o Arabidopsis_thaliana|A1,Vitis_vinifera|V1,Solanum_lycopersicum|L1
INFO: Verified: 7 internal nodes labelled in fam.asr.treefile and present in fam.asr.state;
      topology matches fam.rooted.plain.treefile

fam.asr.treefile:
((((Nepenthes_gracilis|N1:0.243,Drosera_capensis|D1:0.195)Node5:0.129,
   (Cephalotus_follicularis|C1:0.337,Sarracenia_purpurea|S1:0.320)Node6:0.189)Node4:0.246,
   (Utricularia_gibba|U1:0.420,Pinguicula_vulgaris|P1:0.272)Node7:0.204)Node3:0.181,
   (Arabidopsis_thaliana|A1:0.327,Vitis_vinifera|V1:0.231)Node1:0.208,
   Solanum_lycopersicum|L1:0.390)Node2;

$ python3 workflow/scripts/detect_convergence.py -a aln.fa -t fam.asr.treefile \
      -s fam.asr.state --species-yaml config/species.yaml -o fixed.tsv --family rnase_t2
INFO:   7 internal nodes reconstructed
INFO:   all 7 internal tree nodes (and all 9 leaf ancestors) resolved in .state
INFO:   9/9 leaves mapped to species
INFO:   carnivory origin 1: 2 leaf/leaves (Drosera_capensis|D1, Nepenthes_gracilis|N1)
INFO:   carnivory origin 2: 1 leaf/leaves (Cephalotus_follicularis|C1)
INFO:   carnivory origin 3: 1 leaf/leaves (Sarracenia_purpurea|S1)
INFO:   carnivory origin 4: 2 leaf/leaves (Pinguicula_vulgaris|P1, Utricularia_gibba|U1)
INFO:   16 candidate convergent sites before FDR correction
```

7/7 nodes and 9/9 leaf ancestors resolve, versus 1/7 and 2/9 before.

### 5.1 Positive control — a planted convergent substitution is recovered

Column 12 of the simulated alignment was fully conserved (`P` in all nine taxa). `W` was planted in
`Nepenthes_gracilis` (origin 1), `Cephalotus_follicularis` (origin 2) and `Utricularia_gibba`
(origin 4) — one leaf per origin, each with an unmodified sibling so the parent node retains the
ancestral state. The whole two-pass chain was re-run from scratch on the modified alignment:

```
planted: position 12, P -> W, origins 1, 2, 4

rows with n_lineages >= 3:
family    aln_position  ancestral_aa  derived_aa  lineages  n_lineages  min_posterior  p_value      q_value_bh  category
rnase_t2  12            P             W           1,2,4     3           0.9901         2.93e-05     0.000498    full
```

Exact match on position, ancestral state, derived state and the set of origins.

### 5.2 Negative controls — every guard rail fires

```
(a) wrong tree: {family}.rooted.treefile passed instead of .asr.treefile
ValueError: fam.rooted.treefile: 1 internal node(s) have no label, e.g.
  ["<node over ['Drosera_capensis|D1', 'Nepenthes_gracilis|N1']>"]. --tree must be the treefile
  written by the IQ-TREE --ancestral run ...

(b) pass-1 tree from an "-alrt + --ancestral" run (IQ-TREE dropped the root label)
ValueError: p1.treefile: 1 internal node(s) have no label, e.g.
  ["<node over ['Arabidopsis_thaliana|A1', 'Cephalotus_follicularis|C1', ...]...>"] ...

(c)+(d) four leaf labels renamed to unparseable contig IDs (origins 2, 3, 4 destroyed)
WARNING: 4 of 9 tree leaves could not be mapped to any species in config/species.yaml and are
  EXCLUDED from convergence detection: ['contig_00042|C1', 'contig_00043|S1', 'contig_00044|U1',
  'contig_00045|P1']
INFO:   carnivory origin 1: 2 leaf/leaves (Drosera_capensis|D1, Nepenthes_gracilis|N1)
ValueError: Family 'rnase_t2': only 1 independent carnivory origin(s) [1] are represented among
  the mapped leaves, but --min-lineages is 2. No convergent site can possibly be reported.
  Unmapped leaves: ['contig_00042|C1', ...]
```

Under the old code, case (c)/(d) produced an empty TSV and exit 0.

### 5.3 Existing test suite

```
$ python3 -m pytest tests -q
13 passed in 4.76s
```

(`pytest` had to be pip-installed on this laptop; earlier audit entries note it was previously
absent.)

### 5.4 Snakemake DAG

```
$ snakemake -n <all four convergence targets>
Building DAG of jobs...
job                         count
------------------------  -------
infer_tree                      4
root_tree                       4
ancestral_reconstruction        4
detect_convergence              4
```

Correct chaining, four families, no `AmbiguousRuleException`, and no spurious
`family=rnase_t2.asr` jobs (the wildcard constraint works).

---

## 6. Related fix folded in

`workflow/rules/ancestral_structure.smk`'s `extract_ancestor` rule consumed the same mismatched
pair (`{family}.state` + `{family}.rooted.treefile`) and is now repointed to `{family}.asr.state` +
`{family}.asr.treefile`. With a properly labelled tree, `extract_ancestor.py`'s
`if not node_name:` branch — which fabricates a name from a **levelorder index over all nodes**,
a numbering scheme IQ-TREE does not use — can no longer trigger. See §8 for the residual issues in
that script that this task did not address.

---

## 7. Not verified — needs real HPC / real data

1. **Snakemake execution of the new rules.** Snakemake's shell layer fails on this Windows laptop
   for *any* rule, including a two-line `echo` smoke test:
   ```
   rule smoke: shell: "mkdir -p results/ ; echo hello > {output}"
   -> (command exited with non-zero exit code)
   ```
   The DAG was validated by dry run and every shell command in the new rules was executed by hand,
   but the rules themselves have never been run *through* Snakemake. **Run
   `snakemake -n` then a real `--use-conda` run of `phase2` on Hellbender before trusting the
   workflow end to end.**
2. **A pre-existing, unrelated blocker prevents ANY full-workflow Snakemake run**, including dry
   runs, and was found while testing this fix:
   ```
   ChildIOException: File/directory is a child to another output:
   (results/sequences, fetch_all_sequences)
   (results/sequences/aspartic_proteases_a1b_homology/combined.fa, combine_family_sequences)
   ```
   `retrieve.smk`'s `fetch_all_sequences` declares the parent directory of another rule's output.
   This is out of scope for T11 and was **not** fixed here, but it must be fixed before the
   workflow can run at all. Recommend a follow-up task.
3. **Real-data behaviour.** Everything above used a 9-taxon / 300-column simulated alignment.
   Unverified at real scale: IQ-TREE runtime for the ASR pass on 100+ sequences; whether real
   `{Species_name}|{accession}` labels containing unexpected characters survive IQ-TREE's Newick
   round-trip (the `|` character does — verified); whether the outgroup is monophyletic in each real
   gene tree (if not, `run_ancestral.py` will raise, by design — decide per family whether to pass
   `--allow-outgroup-nonmonophyly`).
4. **Model-string round-trip on real ModelFinder output.** `_extract_best_model()` was verified
   against a real report line (`Best-fit model according to BIC: LG`) and complex strings were
   verified to be accepted by `-m` (`LG+F+I+G4` ran fine with `-te --ancestral`), but no real
   family has yet produced a mixture/partition-style verdict. If ModelFinder ever emits something
   `-m` cannot re-parse, the ASR rule will fail loudly — which is the correct behaviour, but budget
   for it.
5. **IQ-TREE binary name on Hellbender.** `config["phylogeny"]["binary"]` defaults to `iqtree3`.
   Confirm the module actually exposes that name.
6. **Statistical model, untouched.** Every candidate site with the same `n_lineages` receives an
   identical p-value, because the Poisson `mu` is site-independent. That is a pre-existing property
   of `_detect_convergent_sites()`, not something this fix introduced or repaired, and it should be
   reviewed before the numbers go into a manuscript.

---

## 8. Known residual issues deliberately left alone

* `extract_ancestor.py` still contains two fabrication fallbacks that should be removed:
  the levelorder-index node-name invention (§6) and a `str.contains(node_name)` substring match on
  the `.state` `Node` column. Repointing the rule makes them unreachable in the normal path, but
  they remain landmines. Recommend a follow-up task.
* `root_tree.py` uses `Bio.Phylo.root_with_outgroup()` without `outgroup_branch_length`. With a
  single outgroup tip this places the root *at* the tip, giving it a zero-length branch
  (`...:0.39030,Solanum_lycopersicum|L1:0.00000):0.00000;` — observed). The ASR pass re-optimises
  branch lengths so the ancestral states are unaffected, but the `{family}.rooted.treefile` used
  for figures is misleading.
* `CLAUDE.md` §3 vs `config/enzyme_families.yaml` tier-list drift was being fixed concurrently by
  another task (T14); not touched here.

---

## 9. Sources

* IQ-TREE Command Reference — <https://iqtree.github.io/doc/Command-Reference>
* IQ-TREE Tutorial — <https://iqtree.github.io/doc/Tutorial>
* IQ-TREE FAQ — <https://iqtree.github.io/doc/Frequently-Asked-Questions>
* IQ-TREE Rootstrap / rooted trees — <https://iqtree.github.io/doc/Rootstrap>
* Google Group, "Question on ancestral sequence reconstruction with IQTREE" (Minh Bui, Georg
  Hochberg) — <https://groups.google.com/g/iqtree/c/jOn0BUN2dGI>
* Cibiv/IQ-TREE issue #9, "Ancestral sequence reconstruction" (bqminh) —
  <https://github.com/Cibiv/IQ-TREE/issues/9>
* IQ-TREE source: `utils/tools.cpp`, `main/treetesting.cpp`, `main/phyloanalysis.cpp`,
  `tree/mtree.cpp` — <https://github.com/iqtree/iqtree2> (default branch `master`)
* IQ-TREE 3 releases — <https://github.com/iqtree/iqtree3/releases/tag/v3.1.3>;
  Windows asset <https://github.com/iqtree/iqtree3/releases/download/v3.1.3/iqtree-3.1.3-Windows.zip>
* `ancseq` (published `-te -asr` usage) — <https://github.com/YuSugihara/ancseq>
* Naser-Khdour S, Minh BQ, Lanfear R. (2021) Assessing Confidence in Root Placement on
  Phylogenies. *Systematic Biology* — <https://doi.org/10.1093/sysbio/syab067>
