# CLAUDE.md — CarnivorEnzyme (Hellbender working copy)

> **This file is tracked in git as of 2026-09-07** (it was previously gitignored, so the laptop and
> the cluster each carried a separate untracked copy). It currently covers *operating on
> Hellbender* only.
>
> **Intended end state: one unified, scoped CLAUDE.md.** The laptop's copy holds the science
> rationale and the Method Justification Table; that content should be **merged into this file**
> rather than kept as a second divergent copy. On the laptop's first pull, merge — do not let it
> overwrite, and do not create a parallel file. See PROJECT_STATUS.md §0.

## Where you are

This checkout lives on Hellbender at
`/cluster/VAST/mendozacozatld-lab/PierceTaylor/CarnivorEnzyme` (the lab's 5 TB VAST share).
Sessions usually run inside an Open OnDemand VS Code job — a **compute node with only ~4 CPUs**,
not a place to run the pipeline directly.

## Hellbender rules

**[docs/HELLBENDER_RULES.md](docs/HELLBENDER_RULES.md) is the authority.** Read it before running
anything. Summary of the parts that bite most often:

- **Never run compute on a login node** — submit it (`sbatch`, or `snakemake --profile config/slurm`).
- **Nothing on this cluster is backed up.** Deleting `results/` is unrecoverable.
- **Maintenance is the 2nd Tuesday**; check `scontrol show reservation` for the exact window, which
  covers the login and OnDemand nodes too.
- **Never pipe `module load`** — the subshell discards the `PATH` change.
- **The cluster's snakemake module cannot submit to SLURM.** Use the project env:
  ```bash
  export PATH=/mnt/pixstor/data/$USER/miniconda/envs/carnivor-smk/bin:$PATH
  export SNAKEMAKE_CONDA_PREFIX=/mnt/pixstor/data/$USER/carnivor-conda
  ```
- **Conda envs must go on pixstor, not the VAST share** — VAST is NFS and rejects `:` in
  filenames, which silently rolls back any env containing perl.

These are enforced by a `PreToolUse` hook (`.claude/hooks/hellbender_guard.py`). If it blocks a
command, the message names the rule and the fix — read it rather than working around it.

## Running the pipeline

Always target an **output file**, never a bare invocation:

```bash
snakemake --profile config/slurm --use-conda -j 20 \
  results/convergence/<family>.convergent_sites.tsv
```

- Bare `snakemake` resolves to `rule phase1` and exits 0 after only downloading sequences
  (known wart, TODO T0.5.5). `all` fails on stub rules.
- `-j` under the SLURM profile is *concurrent cluster jobs*, not local cores.
- Run the driver itself under `sbatch` for anything long — an OnDemand session dies at its walltime
  and takes the Snakemake process with it.
- After an `rsync` of `results/`, run `--touch` first or mtimes will trigger a full NCBI re-fetch.

## Project conventions

- **Never add Claude as a git co-author.** No `Co-Authored-By:` trailer and no session links, in any
  commit or PR, on any repo. Hard rule.
- Audit docs are numbered (`audit/NN_*.md`) and must carry an explicit **"NOT verified"** section
  alongside what was verified by execution. Add a matching `audit/CHANGELOG.md` entry.
- Tier gating: don't start Phase 3 (structure prediction), FoldX, or electrostatics until
  [TODO.md](TODO.md) Tier 0 and Tier 1 close.
- **Exit code 0 is not a pass criterion** for the phylogeny→convergence chain. The bug class in
  `audit/11_ancestral_reconstruction_architecture_fix.md` fails *silently*, producing empty or
  fabricated output. Validate contents against the alignment by hand.
