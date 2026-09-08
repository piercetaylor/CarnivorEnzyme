#!/bin/bash
#SBATCH --job-name=carnivor-smk
#SBATCH --partition=general
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/slurm/driver-%j.out
#SBATCH --error=logs/slurm/driver-%j.err
#
# Snakemake driver job. This process only builds the DAG and submits/waits on the
# real work, so it needs almost no resources itself -- but it must outlive the
# interactive session that launched it. An Open OnDemand session dies at its
# walltime and takes any Snakemake process running inside it down too.
#
# Usage:  sbatch workflow/scripts/run_snakemake.sbatch <target> [extra snakemake args]
#
# Always pass an explicit output file as <target>. A bare `snakemake` resolves to
# `rule phase1` and exits 0 after only downloading sequences (TODO T0.5.5); `all`
# fails on the stub rules.

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "error: no target given. Pass an output file, e.g." >&2
  echo "  sbatch $0 results/convergence/chitinases_gh19_class_iv.convergent_sites.tsv" >&2
  exit 2
fi

# The cluster's snakemake module has no SLURM executor plugin and the env is
# read-only, so use the project env (see docs/HELLBENDER_RULES.md).
SMK_ENV="${SMK_ENV:-/mnt/pixstor/data/$USER/miniconda/envs/carnivor-smk}"
if [[ ! -x "$SMK_ENV/bin/snakemake" ]]; then
  echo "error: no snakemake at $SMK_ENV/bin/snakemake" >&2
  echo "build it with:  mamba create -p $SMK_ENV -c conda-forge -c bioconda \\" >&2
  echo "                  snakemake=8 snakemake-executor-plugin-slurm conda mamba" >&2
  exit 3
fi
export PATH="$SMK_ENV/bin:$PATH"

# Snakemake's per-rule conda envs MUST NOT live under .snakemake/ on this repo's
# filesystem: the lab share is NFS and rejects ':' in filenames, so installing
# perl (a dependency of orthofinder/trimal) dies with
#   [Errno 22] Invalid argument: '.../man/man3/App::Cpan.3'
# and conda rolls the whole environment back. pixstor is GPFS and accepts them.
# This also keeps ~2.5 GB of envs off the 5 TB group allocation.
export SNAKEMAKE_CONDA_PREFIX="${SNAKEMAKE_CONDA_PREFIX:-/mnt/pixstor/data/$USER/carnivor-conda}"
mkdir -p "$SNAKEMAKE_CONDA_PREFIX"

echo "host        : $(hostname)"
echo "driver job  : ${SLURM_JOB_ID:-<none>}"
echo "snakemake   : $(command -v snakemake) ($(snakemake --version))"
echo "conda-prefix: $SNAKEMAKE_CONDA_PREFIX"
echo "target      : $1"
echo "started     : $(date '+%F %T %Z')"
echo

snakemake --profile config/slurm --use-conda "$@"
rc=$?

echo
echo "finished    : $(date '+%F %T %Z')  (exit $rc)"
echo "NOTE: exit 0 is NOT a pass criterion for the convergence chain -- the bug"
echo "      class in audit/11_*.md fails silently. Validate the output contents."
exit $rc
