# Hellbender operating rules

Rules for running CarnivorEnzyme on **Hellbender** (University of Missouri, ITRSS). These are
cluster policy, not project preference — they apply to every user and every clone. Sources are
listed at the bottom.

`CLAUDE.md` and `.claude/` are gitignored in this repo, so this tracked file is the shareable copy;
the local `CLAUDE.md` should point here rather than duplicate it.

## The rules

| # | Rule | Source |
|---|------|--------|
| **R1** | **No compute on login nodes.** "Under no circumstances should your code be running on the login node." Submit with `sbatch`, or take an allocation with `salloc`. | ITRSS wiki |
| **R2** | **Partition walltime caps.** `general` / `requeue` / `gpu` / `logical_cpu` ≤ 2 days; `interactive` ≤ 4 h (and does **not** accept batch jobs); `rss-class` ≤ 8 h; investor `priority` ≤ 28 days. | `sinfo`, ITRSS wiki |
| **R3** | **No account sharing.** "Each user must use their own account to access RSS resources. Account sharing is prohibited." There is no root/sudo. | ITRSS wiki |
| **R4** | **Storage quotas.** `/home/$USER` 50 GB; `/home/$USER/data` 500 GB; this lab's VAST group share is 5 TB. | ITRSS wiki |
| **R5** | **Nothing is backed up.** "None of the cluster attached storage available to users is backed up in any way by us." Deletion is irreversible — there are no snapshots. | ITRSS wiki |
| **R6** | **Data Class 1 and 2 only.** Cluster storage is "limited to Data Class 1 and 2 as defined by UM System DCL". DCL3/DCL4 (PHI, HIPAA, identifiable human-subject data) must not be stored here. | ITRSS wiki |
| **R7** | **Maintenance window.** The 2nd Tuesday of each month. Jobs that cannot finish before it starts are held until it ends. | ITRSS wiki |
| **R8** | **Licensed software needs approval.** Procurement goes through ITSRQ. Relevant here: FoldX 5.1 is licensed. | ITRSS wiki |

R7's exact window is published as a SLURM reservation and is authoritative — prefer it over
counting Tuesdays:

```bash
scontrol show reservation        # e.g. Sept_Maint StartTime=2026-09-08T08:00:00 ... Flags=MAINT
```

Note the reservation covers **all nodes**, including the login and Open OnDemand hosts — an
interactive session is killed too, not just batch jobs.

## How these are enforced

`.claude/hooks/hellbender_guard.py` runs as a `PreToolUse` hook on every Bash command and returns
`deny` for hard rules, `ask` for anything risking irreversible data loss or resting on a heuristic.
`.claude/hooks/session_context.sh` prints placement and quota context at session start. Behaviour
tests:

```bash
python3 .claude/hooks/test_hellbender_guard.py     # 24 cases
```

The guard **fails open** — an unexpected exception is reported on stderr and the command proceeds,
because a crashing hook that blocks every Bash call would be worse than a missed check. It is a
guardrail against mistakes, not a security boundary against a determined actor.

## Cluster gotchas (learned the hard way)

**Never pipe `module load`.** `module load x 2>&1 | tail -3` runs the function in a *subshell*, so
the `PATH` change is discarded and the tool appears missing. In non-interactive shells also source
Lmod first:

```bash
source /usr/share/lmod/lmod/init/bash
module load mamba/v2.4.0          # NOT: module load mamba/v2.4.0 | tail
```

**The cluster Snakemake module cannot submit to SLURM.** `module load snakemake/v8.27.1` gives
`--executor {local,dryrun,touch}` — `snakemake-executor-plugin-slurm` is absent and the env is
read-only, so `--profile config/slurm` (which declares `executor: slurm`) fails at plugin load.
Use the project env instead, which also carries a `conda` frontend for `--use-conda`:

```bash
export PATH=/mnt/pixstor/data/$USER/miniconda/envs/carnivor-smk/bin:$PATH
snakemake --profile config/slurm --use-conda <target>
```

Built with:

```bash
mamba create -p /mnt/pixstor/data/$USER/miniconda/envs/carnivor-smk \
  -c conda-forge -c bioconda snakemake=8 snakemake-executor-plugin-slurm conda mamba
```

**System `python3` is 3.6.8.** `subprocess.run(capture_output=..., text=...)` is 3.7+ and raises
`TypeError` here. Hook code must stay 3.6-compatible (`stdout=subprocess.PIPE`,
`universal_newlines=True`), and must not hide the failure behind a bare `except: pass`.

**Quotas are not queryable.** `quota` reports "none"; there is no `mmlsquota` or `lfs`. `df` on
`$HOME` shows the whole shared 60 TB GPFS mount, not the 50 GB personal quota. `du -sh $HOME` is
the only accurate check and takes minutes. Keep conda off `$HOME` — `~/.condarc` already redirects
`envs_dirs`/`pkgs_dirs` to `/mnt/pixstor/data/$USER/miniconda`.

**The lab VAST share is NFS and rejects `:` in filenames.** Installing perl (pulled in by
orthofinder/trimal) fails there with `[Errno 22] Invalid argument: '.../man/man3/App::Cpan.3'`, and
conda rolls the entire environment back. Snakemake's per-rule envs must therefore live on pixstor
(GPFS), not under `.snakemake/` in the repo:

```bash
export SNAKEMAKE_CONDA_PREFIX=/mnt/pixstor/data/$USER/carnivor-conda
```

`workflow/scripts/run_snakemake.sbatch` sets this automatically. Watch for it: conda reports the
rollback but the surrounding shell can still exit 0, especially if the output is piped (`| tail`
discards the exit code unless `set -o pipefail` is active).

**Keep large writes off `$HOME`.** Conda envs, Snakemake's `.snakemake/conda`, and results belong on
pixstor or the lab VAST share.

## Sources

- [Hellbender documentation — ITRSS Wiki](https://docs.itrss.umsystem.edu/pub/hpc/hellbender)
- [Getting Started with HPC — ITRSS Wiki](https://docs.itrss.umsystem.edu/pub/hpc/start)
- [Hellbender HPC Cluster Information and Policies — DoIT KB 1310](https://tdx.umsystem.edu/TDClient/36/DoIT/KB/ArticleDet?ID=1310)
- [UM System Collected Rules 110.005 — Acceptable Use Policy](https://www.umsystem.edu/ums/rules/collected_rules/facilities/ch110/110.005_acceptable_use_policy)
