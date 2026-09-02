# 16 — Snakemake DAG: ChildIOException fix (TODO T0.2)

> Date: 2026-09-02
> Task: TODO.md T0.2 — "Fix the `ChildIOException` Snakemake DAG bug"
> Blocks: TODO T0.1 (the Hellbender phylogeny→convergence verification run), and with it
> everything downstream of Phase 2.
> Files changed: `workflow/rules/retrieve.smk`, `workflow/rules/alignment.smk`,
> `workflow/rules/esm2.smk`, `workflow/scripts/score_esm2.py`

---

## 1. The bug, as reproduced

Running the exact command TODO T0.1 calls for, before any change:

```
$ python -m snakemake -n results/convergence/chitinases_gh19_class_iv.convergent_sites.tsv
Building DAG of jobs...
ChildIOException:
File/directory is a child to another output:
(WindowsPath('.../results/sequences'), fetch_all_sequences)
(WindowsPath('.../results/sequences/chitinases_gh19_class_iv/combined.fa'), combine_family_sequences)
```

This fires during DAG construction, before any job is scheduled, so **no dry run past Phase 1 was
possible** — including the Tier 0 verification run the whole project is currently gated on.

## 2. Root cause

`retrieve.smk` declared two outputs in a parent/child relationship:

| Rule | Output |
|---|---|
| `fetch_all_sequences` (checkpoint) | `directory("results/sequences/")` |
| `combine_family_sequences` | `results/sequences/{family}/combined.fa` |

Snakemake forbids this: if one rule declares a directory as its output, it owns that subtree, and
no other rule may write into it. The directory output is correct and necessary here — the checkpoint
exists precisely because the set of downloaded FASTAs is not knowable until runtime. The error was
placing a *derived* artefact (the per-family concatenation) inside the *downloaded* artefact's tree.

## 3. The fix

Move the derived file out of the checkpoint's subtree. `combined.fa` now lives at
**`results/family_fasta/{family}.combined.fa`**.

- `retrieve.smk` — `combine_family_sequences.output` repointed; module docstring gained a path note
  explaining why nothing may be written under `results/sequences/`.
- `alignment.smk` — `align_family.input` repointed.
- `_get_species_fastas()` — kept the filename exclusion as a **defensive guard**, widened to cover
  `all_sequences.fa` too. It is no longer load-bearing (the file it excluded now lives elsewhere),
  but a `combined.fa` left over from a pre-fix run would otherwise be globbed up and concatenated
  into its own successor. There are no such leftovers in the current working tree — checked.

## 4. Second finding, fixed: `esm2.smk` had an input no rule produced

While verifying the DAG, `phase3a` turned out to be unbuildable for an unrelated reason.
`esm2.smk` declared:

```python
sequences="results/sequences/{family}/all_sequences.fa",
```

**No rule anywhere in the workflow produces `all_sequences.fa`.** `combine_family_sequences`
produces `combined.fa`. The two names refer to the same thing — the per-family concatenation — and
the mismatch meant Phase 3A would have failed with `MissingInputException` the moment the
ChildIOException stopped masking it.

Repointed to the same `results/family_fasta/{family}.combined.fa`, along with the stale path in
`score_esm2.py`'s docstring. This is a scope expansion beyond T0.2 as written; it is included
because it is the same class of defect in the same dependency chain, and because T0.2's own pass
criterion ("dry run reaches DAG construction") could not be honestly claimed for `phase3a` while it
stood.

## 5. Third finding — FLAGGED, NOT FIXED: the default target is `phase1`, not `all`

`snakemake -n` with no target builds **2 jobs and stops at Phase 1**:

```
job                    count
fetch_all_sequences        1
phase1                     1
total                      2
```

Snakemake takes the first rule in the file as the default target. `rule phase1` is defined at
Snakefile:83, `rule all` at Snakefile:141. So the bare invocation resolves to `phase1`.

The Snakefile's own usage docstring says otherwise:

```
  # Run full pipeline (all families)
  snakemake --use-conda --cores 8
```

That command does not run the full pipeline; it downloads sequences and exits 0. Anyone trusting
the docstring would read a successful Phase-1-only run as a successful full-pipeline run.

**Not fixed here.** The repair is to move `rule all` above the phase targets, or add
`default_target: True` to it — but that changes the default behaviour of the entire workflow, which
is a decision to make deliberately rather than fold into a DAG bugfix. Filed for TODO Tier 0.5.

## 6. Verification

All by execution, not inspection. Snakemake 9.25.2, Python 3.13.14, Windows.

| Target | Before fix | After fix |
|---|---|---|
| `results/convergence/chitinases_gh19_class_iv.convergent_sites.tsv` (the T0.1 target) | `ChildIOException` | **DAG OK, 8 jobs** — `fetch_all_sequences` → `combine_family_sequences` → `align_family` → `trim_alignment` → `infer_tree` → `root_tree` → `ancestral_reconstruction` → `detect_convergence` |
| `phase1` | `ChildIOException` | DAG OK, 2 jobs |
| `phase2` | `ChildIOException` | DAG OK, 26 jobs (4 tier-1 families) |
| `phase3a` | `ChildIOException`, and `MissingInputException` behind it | DAG OK, 34 jobs |
| `phase4a` | `ChildIOException` | `MissingInputException` on `results/structures/` — **expected**, `predict_structure.smk` is a stub (Tier 2) |
| `phase5d` | `ChildIOException` | `MissingInputException` — same cause, expected |
| `all` | `ChildIOException` | `MissingInputException` on `results/atlas/*` — expected, `integrate.smk` is a stub |
| `pytest tests/` | 13 passed | 13 passed |

The chain in row 1 is exactly the four-rule Phase-2 sequence that `audit/11` established and
depends on, reached in the correct order.

## 7. Not verified

- **Nothing was run on Hellbender.** These are dry runs on Windows. Snakemake's shell layer fails
  locally for even a trivial rule (noted in `audit/11`), so actual execution of the chain remains
  TODO T0.1's job.
- **No job was executed at all** — `-n` only builds the DAG. That the DAG is well-formed says
  nothing about whether `fetch_sequences.py`, MAFFT, or IQ-TREE succeed on real inputs.
- **One thing to expect on the cluster:** the dry run schedules `fetch_all_sequences` to re-run,
  with reason "Updated input files: config/enzyme_families.yaml". After rsyncing
  `results/sequences/` up, Snakemake will still want to re-download unless the timestamps are
  reconciled — `snakemake --touch <target>` before the real run, or let it re-fetch if the compute
  node has NCBI access. This is mtime bookkeeping, not a defect.
