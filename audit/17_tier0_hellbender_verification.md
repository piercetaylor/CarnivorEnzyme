# 17 — Tier 0 Hellbender verification run (TODO T0.1)

> Date: 2026-09-07
> Task: TODO.md T0.1 — "Run the real phylogeny→convergence chain on Hellbender"
> Result: **FAIL** — the chain halts at `ancestral_reconstruction`. The target
> `results/convergence/chitinases_gh19_class_iv.convergent_sites.tsv` was **not produced**.
> Files changed: `workflow/scripts/run_snakemake.sbatch` (new),
> `workflow/scripts/validate_convergence_output.py` (new), `docs/HELLBENDER_RULES.md` (new),
> `CLAUDE.md` (new), `.claude/` hooks (new). **No pipeline logic was changed.**

---

## 1. What was run

First execution of the phylogeny→convergence chain on Hellbender against real project data. Repo at
`d6c73a1`, family `chitinases_gh19_class_iv`, sequences staged by rsync (50 FASTAs across 7
families; 9 species files for this family).

```
sbatch workflow/scripts/run_snakemake.sbatch \
  results/convergence/chitinases_gh19_class_iv.convergent_sites.tsv -j 20
```

Driver SLURM job `17171597` on `c001`, 2026-09-07 20:55:13 → 21:01:22 CDT, elapsed **00:07:08**,
final state FAILED (exit 1:0).

## 2. Verified by execution

**IQ-TREE binary — the explicit checklist item.** `config/config.yaml`'s `phylogeny.binary: iqtree3`
is **correct and needed no change**. The bioconda build installs `iqtree3` as the real binary with
`iqtree` as a symlink to it:

```
$ iqtree3 --version
IQ-TREE version 3.1.3 for Linux x86 64-bit built Jul 26 2026
```

Same major.minor.patch as the local Windows 3.1.3 binary the fix was originally developed against.

**Five of seven rules completed on real cluster software**, in the correct order:

| # | Rule | SLURM job | Outcome |
|---|------|-----------|---------|
| 1 | `combine_family_sequences` | 17171651 | OK — 18 sequences from 9 species |
| 2 | `align_family` (MAFFT 7.525) | 17171682 | OK |
| 3 | `trim_alignment` (trimAl 1.4.rev15) | 17171682 | OK — 18 × 268 columns |
| 4 | `infer_tree` (IQ-TREE 3.1.3, WAG+G4, 1000 UFBoot) | 17171735 | OK |
| 5 | `root_tree` | 17171789 | OK |
| 6 | `ancestral_reconstruction` | 17171870 | **FAILED** |
| 7 | `detect_convergence` | — | never ran |

**The DAG is 7 jobs, not the 8 the T0.1 prompt predicted.** `fetch_all_sequences` is correctly
absent: the rsynced FASTAs are current, so Phase-1 retrieval is skipped and no NCBI access is
needed. The `snakemake --touch` step the prompt prescribes as insurance against an mtime-triggered
re-fetch is therefore **unnecessary** — the dry run scheduled no fetch. (audit/16's closing note
predicted the opposite; that prediction did not hold here.)

**No `ChildIOException`** at any point — the T0.2 fix holds on real data.

**Test suite: 13/13 pass** in the built `bioinfo` env.

## 3. Why it failed

`ancestral_reconstruction` aborted with a deliberate, hardened error:

```
ValueError: IQ-TREE reported 'Branch separating outgroup is not found': the outgroup tips
['Arabidopsis_thaliana|BAD44251.1', 'Solanum_lycopersicum|XP_004237833.1',
 'Arabidopsis_thaliana|O24603.1', 'Secale_cereale|Q9FRV0']
do not form a clade in ...rooted.plain.treefile
```

**This is the guard working, not a bug.** The script refused to write an ASR tree it could not
orient, rather than emitting a silently mis-rooted one.

The underlying cause is **not** rooting — it is **paralogy**. The four designated outgroup tips are
scattered across the ML tree, interleaved with carnivore sequences:

- `Solanum_lycopersicum|XP_004237833.1` is **sister to `Sarracenia_purpurea|KAL6971754.1`**
- `Secale_cereale|Q9FRV0` is **sister to the two `Dionaea_muscipula` sequences**
- `Arabidopsis_thaliana|O24603.1` sits inside a clade of Dionaea + Secale + Cephalotus
- only `Arabidopsis_thaliana|BAD44251.1` is basal

Supporting evidence that this is a gene-family scope problem rather than a tree-inference artifact:
the alignment holds **18 sequences from 9 species** — `Nepenthes_gracilis` contributes 4,
`Cephalotus_follicularis` 3, and four more species 2 apiece. Several are paralogs that duplicated
before the species split, so no species-level outgroup can be monophyletic. One branch,
`Cephalotus_follicularis|GAV70878.1` at length **1.73**, is an order of magnitude longer than its
neighbours and is a candidate mis-assignment.

Convergence detection over a tree mixing paralogs would compare non-orthologous sites, which is a
mechanism for producing *exactly* the spurious "convergent substitutions" this project exists to
find. That makes this a scope question, not a flag to switch off.

## 4. The escape hatch, and why it was not used

`run_ancestral.py:203` offers `--allow-outgroup-nonmonophyly`. The error text argues the risk is
narrow — "Ancestral states themselves are unaffected (marginal ASR under a reversible model is
root-independent), but any downstream MRCA lookup on this tree is unreliable" — and that argument
checks out for this consumer specifically: `detect_convergence.py` resolves each leaf's **direct
parent** via `_leaf_parent_names()` (`leaf.up`, lines 242–258) and never calls
`get_common_ancestor`. So the flag would very likely produce a TSV.

It was **not** used, for two reasons:

1. The T0.1 prompt is explicit — do not silently work around a mismatch; report it.
2. Passing the flag treats the symptom. The rooting failure is the detector for a paralogy problem
   that would still be present, unflagged, in every downstream number.

This is a scientific call for Pierce, not a mechanical one. See §6.

## 5. NOT verified

- **`detect_convergence.py` never executed.** The entire premise of T0.1 — that the `9b1a8c5`
  node-matching fix works against real alignments and a real IQ-TREE `.state` file — remains
  **untested**. Nothing in this run touched that code path.
- No `.asr.treefile` or `.asr.state` was produced, so pass criteria 1–4 are all **unevaluated**,
  not passed and not failed.
- The other three Tier-1 families (`purple_acid_phosphatase`, `rnase_t2`,
  `aspartic_proteases_a1b_homology`) were not run, per the prompt's gate. Their outgroup sets have
  **not** been checked for the same non-monophyly — this is likely a family-wide issue, not one
  specific to chitinases.
- `piqtree` is listed in `workflow/envs/bioinfo.yaml` with the comment "IQ-TREE 3 Python bindings
  for tree manipulation in detect_convergence.py". Nothing in `workflow/` or `tests/` imports it;
  `detect_convergence.py` uses `ete3`. The comment is stale (the package installs fine).
- The T0.1 prompt names the output column `alignment_position`; the code writes **`aln_position`**
  (`detect_convergence.py:536`). Doc drift only.

## 6. Agreed direction — recover as many common ancestors as possible

**Decision (Pierce, 2026-09-07): do not re-run yet.** The next attempt is re-framed away from
"make the outgroup monophyletic" toward **maximising the number of internal nodes (common
ancestors) that can be reconstructed and trusted for this gene tree.**

The rationale is that outgroup monophyly is not actually what this analysis needs:

- **Marginal ASR is root-independent under a reversible model.** The run used WAG+G4, which is
  reversible, so per-node marginal ancestral states do not depend on where the tree is rooted. The
  error message in `run_ancestral.py` says as much.
- **`detect_convergence.py` never performs an MRCA lookup.** It resolves each leaf's *direct
  parent* through `_leaf_parent_names()` (`leaf.up`, lines 242–258) and requires only that every
  such parent appear in the `.state` file (the check at line 284). It never calls
  `get_common_ancestor`.

So the binding requirement is **per-node state coverage**, not orientation. Outgroup monophyly is a
convenience for interpreting *direction* of change, not a prerequisite for reconstructing states.

### Planned work

1. **Quantify the recoverable ancestor set.** Re-run `ancestral_reconstruction` with
   `--allow-outgroup-nonmonophyly` and report, as a table: total internal nodes in
   `.asr.treefile`; how many appear in `.asr.state`; how many are the direct parent of at least one
   leaf; and how many of those parents subtend leaves from ≥2 distinct `carnivory_origin` values
   (the minimum for a convergence call at `min_lineages: 2`). For 18 tips the ceiling is 16–17
   internal nodes; the useful figure is how many are *leaf-parents with origin diversity*.
2. **Record per-node posterior support.** `.state` carries per-site posteriors; summarise the
   distribution per node so low-confidence ancestors can be excluded rather than silently averaged
   in. `config.yaml convergence.posterior_threshold` is already 0.8.
3. **Separate orientation from reconstruction in the code.** `run_ancestral.py` currently makes
   non-monophyly a hard failure for *all* consumers. Given point 2 above, the better shape is to
   reconstruct states regardless, and mark the output tree as "not outgroup-oriented" so that any
   future consumer needing MRCA/polarity fails loudly while `detect_convergence` proceeds. That is
   a code change and belongs in its own task, not a flag flipped at the command line.
4. **Report which ancestors are *not* recoverable and why** — nodes whose descendants are all one
   origin, or whose posteriors fall below threshold. That set is the honest limit of what this gene
   tree can support, and it belongs in the eventual results.

### Still true, and still a separate problem

The paralogy remains: 18 sequences from 9 species, with `Solanum_lycopersicum` sister to
`Sarracenia_purpurea`, `Secale_cereale` sister to the two `Dionaea_muscipula` sequences, and
`Cephalotus_follicularis|GAV70878.1` on a 1.73 branch. Maximising ancestor recovery does **not**
resolve that; a substitution shared between paralogous lineages is not evidence of convergent
adaptation. Family re-scoping to a true orthogroup — anticipated by
`archive/ORTHOLOGY_AND_FAMILY_SCOPE_2026-05-12.md` — stays on the table as separate work, and any
convergence numbers produced before it lands must be labelled provisional.

## 7. Infrastructure fixed along the way

Three Hellbender-specific blockers stood between the repo and a first run. None are in the T0.1
prompt. All are now documented in `docs/HELLBENDER_RULES.md` and enforced where possible.

1. **The cluster's Snakemake module cannot submit to SLURM.** `module load snakemake/v8.27.1`
   exposes `--executor {local,dryrun,touch}`; `snakemake-executor-plugin-slurm` is absent and the
   env is read-only, so `--profile config/slurm` (which declares `executor: slurm`) fails at plugin
   load. Fixed by building a project driver env on pixstor:
   `mamba create -p /mnt/pixstor/data/$USER/miniconda/envs/carnivor-smk -c conda-forge -c bioconda
   snakemake=8 snakemake-executor-plugin-slurm conda mamba` (Snakemake 8.30.0).
2. **The lab VAST share is NFS and rejects `:` in filenames.** Building `bioinfo.yaml` under
   `.snakemake/` died installing perl —
   `[Errno 22] Invalid argument: '.../man/man3/App::Cpan.3'` — and conda rolled the **entire**
   environment back. The surrounding shell still reported exit 0 because the output was piped
   (`| tail` discards the status without `set -o pipefail`). Fixed by putting per-rule envs on
   pixstor (GPFS) via `SNAKEMAKE_CONDA_PREFIX=/mnt/pixstor/data/$USER/carnivor-conda`.
   `workflow/scripts/run_snakemake.sbatch` sets this; the PreToolUse hook now denies a
   `--use-conda` run on NFS without a prefix.
3. **`module load` piped into another command silently no-ops** — the pipeline runs the shell
   function in a subshell, so the `PATH` change is discarded and the tool appears missing.
   Non-interactive shells also need `source /usr/share/lmod/lmod/init/bash` first.

Non-blocking notes: Snakemake logs "Unable to guess SLURM account. Trying to proceed without" —
harmless here, every child job submitted and ran. Neither env ships `pytest`; it was pip-installed
into the `bioinfo` env to run the suite.

## 8. Reproducing this run

```bash
source /usr/share/lmod/lmod/init/bash          # do NOT pipe module commands
export PATH=/mnt/pixstor/data/$USER/miniconda/envs/carnivor-smk/bin:$PATH
export SNAKEMAKE_CONDA_PREFIX=/mnt/pixstor/data/$USER/carnivor-conda

snakemake --profile config/slurm --use-conda -n \
  results/convergence/chitinases_gh19_class_iv.convergent_sites.tsv     # 7 jobs

sbatch workflow/scripts/run_snakemake.sbatch \
  results/convergence/chitinases_gh19_class_iv.convergent_sites.tsv -j 20
```

Once a TSV exists, validate it rather than trusting the exit code:

```bash
python3 workflow/scripts/validate_convergence_output.py chitinases_gh19_class_iv --rows 3
```

That script re-derives each sampled row from the trimmed alignment and the `.state` file, and
checks that every internal tree label resolves — the four T0.1 pass criteria, mechanised.
