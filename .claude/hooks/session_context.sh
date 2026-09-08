#!/usr/bin/env bash
# SessionStart hook: report where this session is running on Hellbender, and how
# close the project is to its storage quotas. Placement context up front means an
# agent never has to guess whether it is on a login node.
#
# Deliberately cheap: `df` not `du`. A recursive `du` of /home takes minutes here
# and would stall every session start.

set -uo pipefail

host="$(hostname)"
job="${SLURM_JOB_ID:-${SLURM_JOBID:-}}"
cpus="${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo '?')}"

echo "── Hellbender session context ───────────────────────────────"
if [[ -n "$job" ]]; then
  echo "  node        : $host  (SLURM job $job, ${cpus} CPU)"
  echo "  placement   : inside an allocation — compute is allowed here"
else
  echo "  node        : $host  (no SLURM allocation)"
  echo "  placement   : LOGIN NODE — do not run compute. Use sbatch/salloc,"
  echo "                or: snakemake --profile config/slurm --use-conda <target>"
fi

vast="$(df -h --output=used,size,pcent /cluster/VAST/mendozacozatld-lab 2>/dev/null | tail -1 | tr -s ' ')"
[[ -n "${vast// /}" ]] && echo "  lab VAST    :${vast} used (5 TB group allocation)"

# NOTE: $HOME sits on the shared GPFS mount, so `df` reports the whole 60 TB
# filesystem, not this user's 50 GB quota -- and Hellbender exposes no `quota`,
# `mmlsquota` or `lfs` command to query the real figure. Reporting df here would
# be actively misleading, so state the policy instead. `du -sh $HOME` is the only
# accurate check and takes minutes; run it manually if a quota error appears.
echo "  \$HOME quota : 50 GB (not queryable here — conda is redirected to pixstor"
echo "                via ~/.condarc, which is what usually fills it)"

if command -v squeue >/dev/null 2>&1; then
  n="$(squeue -h -u "$USER" 2>/dev/null | wc -l | tr -d ' ')"
  echo "  your jobs   : $n queued/running"
fi

# R7: report the exact maintenance outage from SLURM's MAINT reservation, which
# is authoritative and survives a rescheduled window; fall back to the documented
# 2nd-Tuesday rule if scontrol is unavailable.
python3 - <<'MAINT_PY' 2>/dev/null
import calendar, datetime, re, subprocess
best = None
try:
    out = subprocess.run(["scontrol", "show", "reservation"],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         universal_newlines=True, timeout=5).stdout
    for block in out.split("ReservationName="):
        if "MAINT" not in block:
            continue
        m = re.search(r"StartTime=(\S+)", block)
        if not m:
            continue
        try:
            st = datetime.datetime.strptime(m.group(1), "%Y-%m-%dT%H:%M:%S")
        except ValueError:
            continue
        if st > datetime.datetime.now() and (best is None or st < best):
            best = st
except Exception:
    pass
if best:
    hrs = (best - datetime.datetime.now()).total_seconds() / 3600
    print(f"  maintenance : {best:%Y-%m-%d %H:%M} — in {hrs:.1f} h; jobs that cannot")
    print("                finish before it starts will be held until it ends")
else:
    today = datetime.date.today()
    for off in (0, 1):
        y = today.year + (today.month - 1 + off) // 12
        m = (today.month - 1 + off) % 12 + 1
        tues = [d for d in calendar.Calendar().itermonthdates(y, m)
                if d.month == m and d.weekday() == 1]
        if len(tues) >= 2 and tues[1] >= today:
            print(f"  maintenance : {tues[1]} (2nd Tuesday) — long jobs will be held")
            break
MAINT_PY

echo "  snakemake   : use /mnt/pixstor/data/\$USER/miniconda/envs/carnivor-smk"
echo "                (the cluster snakemake module has no SLURM executor plugin)"
echo "─────────────────────────────────────────────────────────────"
