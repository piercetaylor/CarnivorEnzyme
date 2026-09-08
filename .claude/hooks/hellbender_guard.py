#!/usr/bin/env python3
"""PreToolUse guard enforcing Hellbender (Mizzou ITRSS) cluster rules on Bash commands.

Rules enforced (see CLAUDE.md "Hellbender operating rules" for the prose version):

  R1  No compute on login nodes -- "Under no circumstances should your code be
      running on the login node."
  R2  Partition walltime caps (read live from `sinfo`, static fallback below).
  R3  No account sharing -- each user uses their own account.
  R5  Cluster storage is NOT backed up; destructive ops need confirmation.
  R6  Cluster storage is limited to UM System Data Class 1 and 2.
  R7  Maintenance window: 2nd Tuesday of each month.

Plus two project-local rules:
  * Snakemake that actually executes must be routed through the SLURM profile.
  * Thread/core requests must not exceed the current allocation's CPU count.

Decisions: hard rules -> "deny"; anything that risks irrecoverable data loss or
is only a heuristic -> "ask" (a human confirms). Clean commands produce no
output, so the harness applies its normal permission logic.

Fails OPEN: an unexpected exception here must not brick every Bash call. The
failure is announced on stderr so it is visible rather than silent.
"""

import json
import os
import re
import shlex
import socket
import subprocess
import sys
from datetime import date, datetime, timedelta

# --- R2: static fallback if `sinfo` is unavailable (minutes) ---------------
STATIC_PARTITION_MAX_MIN = {
    "general": 2 * 24 * 60,
    "requeue": 2 * 24 * 60,
    "gpu": 2 * 24 * 60,
    "gpu_requeue": 2 * 24 * 60,
    "logical_cpu": 2 * 24 * 60,
    "logical_cpu2": 2 * 24 * 60,
    "interactive": 4 * 60,
    "rss-class": 8 * 60,
}

# --- R1: things that must never run outside a job allocation ---------------
HEAVY_BINARIES = {
    "iqtree", "iqtree2", "iqtree3", "raxml", "raxml-ng", "raxmlHPC",
    "mafft", "mafft-linsi", "mafft-einsi", "muscle", "trimal", "clustalo",
    "hmmsearch", "hmmscan", "hmmbuild", "jackhmmer", "phmmer",
    "orthofinder", "orthofinder.py", "diamond", "mmseqs", "cd-hit",
    "blastp", "blastn", "blastx", "tblastn", "psiblast", "makeblastdb",
    "gmx", "gmx_mpi", "foldx", "FoldX", "vina", "smina", "salmon",
    "apbs", "pdb2pqr", "orca", "colabfold_batch", "alphafold",
}
VERSION_FLAGS = {"--version", "-version", "--help", "-h", "-V", "--v"}

# Thread-count flags, per program, that consume local CPUs.
THREAD_FLAGS = {
    "orthofinder": {"-t", "-a"}, "orthofinder.py": {"-t", "-a"},
    "iqtree": {"-T", "-nt"}, "iqtree2": {"-T", "-nt"}, "iqtree3": {"-T", "-nt"},
    "mafft": {"--thread"}, "mafft-linsi": {"--thread"},
    "hmmsearch": {"--cpu"}, "hmmscan": {"--cpu"}, "jackhmmer": {"--cpu"},
    "diamond": {"-p", "--threads"}, "mmseqs": {"--threads"},
    "salmon": {"-p", "--threads"}, "trimal": set(),
    "blastp": {"-num_threads"}, "blastn": {"-num_threads"},
    "blastx": {"-num_threads"}, "tblastn": {"-num_threads"},
}

# Logins that are a service convention rather than a shared human account, so R3
# does not apply: `git@github.com` authenticates with the caller's own SSH key.
SERVICE_LOGINS = {"git", "hg", "svn", "aur"}
FORGE_HOSTS = {
    "github.com", "gitlab.com", "bitbucket.org", "codeberg.org",
    "git.sr.ht", "ssh.github.com", "altssh.gitlab.com",
}

# --- R5: storage with no backup; destructive ops here need a human ---------
PROTECTED_PREFIXES = [
    "/cluster/VAST/mendozacozatld-lab",
    "/mnt/pixstor/data/",
    os.path.expanduser("~/data"),
]
PROTECTED_RELATIVE = ("results/", "results", "logs/", "resources/")

# --- R6: Data Class 3/4 tripwire (heuristic, not real DLP) -----------------
DCL_PATTERN = re.compile(r"(?i)(^|[_\-./])(phi|hipaa|ssn|mrn|pii|hcpcs)([_\-.]|$)")

SNAKEMAKE_NO_EXEC = {
    "-n", "--dry-run", "--dryrun", "--dag", "--rulegraph", "--filegraph",
    "--summary", "--detailed-summary", "--list", "--list-target-rules",
    "--list-code-changes", "--list-input-changes", "--list-params-changes",
    "--unlock", "--cleanup-metadata", "--conda-create-envs-only",
    "--containerize", "--report", "--generate-unit-tests", "--version", "--help",
}

CONTROL_TOKENS = {"&&", "||", ";", "|", "&", "\n"}


def decide(decision, reason):
    """Emit a PreToolUse decision and exit."""
    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": decision,
            "permissionDecisionReason": reason,
        }
    }))
    sys.exit(0)


def partition_limits():
    """Live partition walltime caps from sinfo, falling back to the static table."""
    limits = dict(STATIC_PARTITION_MAX_MIN)
    try:
        out = subprocess.run(
            ["sinfo", "-h", "-o", "%P|%l"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, timeout=5,
        ).stdout
        for line in out.splitlines():
            if "|" not in line:
                continue
            name, tl = line.split("|", 1)
            name = name.strip().rstrip("*")
            mins = parse_slurm_time(tl.strip())
            if name and mins:
                limits[name] = mins
    except Exception as exc:
        print("[hellbender_guard] sinfo lookup failed (%r); using static "
              "partition table" % (exc,), file=sys.stderr)
    return limits


def parse_slurm_time(text):
    """SLURM walltime -> minutes. Handles m, m:s, h:m:s, d-h, d-h:m, d-h:m:s."""
    if not text or text in ("UNLIMITED", "INFINITE", "infinite", "NOT_SET"):
        return None
    text = text.strip()
    days = 0
    if "-" in text:
        d, _, text = text.partition("-")
        try:
            days = int(d)
        except ValueError:
            return None
    parts = text.split(":") if text else ["0"]
    try:
        nums = [int(p) for p in parts]
    except ValueError:
        return None
    if len(nums) == 1:
        h, m, s = 0, nums[0], 0
    elif len(nums) == 2:
        h, m, s = 0, nums[0], nums[1]
    elif len(nums) == 3:
        h, m, s = nums
    else:
        return None
    return days * 1440 + h * 60 + m + (1 if s else 0)


def maintenance_start():
    """Exact start of the next maintenance outage.

    SLURM carries a MAINT reservation (e.g. `Sept_Maint`, StartTime=
    2026-09-08T08:00:00, Flags=MAINT,...,ALL_NODES), which is authoritative and
    survives a rescheduled window. Fall back to the documented 2nd-Tuesday rule
    only when scontrol is unavailable.

    Returns a datetime, or None.
    """
    try:
        out = subprocess.run(
            ["scontrol", "show", "reservation"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, timeout=5,
        ).stdout
        best = None
        for block in out.split("ReservationName="):
            if "MAINT" not in block:
                continue
            m = re.search(r"StartTime=(\S+)", block)
            if not m:
                continue
            try:
                st = datetime.strptime(m.group(1), "%Y-%m-%dT%H:%M:%S")
            except ValueError:
                continue
            if st > datetime.now() and (best is None or st < best):
                best = st
        if best:
            return best
    except Exception as exc:
        print("[hellbender_guard] scontrol lookup failed (%r); falling back to "
              "the 2nd-Tuesday rule" % (exc,), file=sys.stderr)
    d = next_maintenance_tuesday()
    return datetime.combine(d, datetime.min.time()) if d else None


def next_maintenance_tuesday(today=None):
    """Documented fallback: maintenance is the 2nd Tuesday of each month (R7)."""
    today = today or date.today()
    for month_offset in (0, 1, 2):
        y = today.year + (today.month - 1 + month_offset) // 12
        m = (today.month - 1 + month_offset) % 12 + 1
        tuesdays = [
            date(y, m, d)
            for d in range(1, 32)
            if _valid(y, m, d) and date(y, m, d).weekday() == 1
        ]
        if len(tuesdays) >= 2 and tuesdays[1] >= today:
            return tuesdays[1]
    return None


def _valid(y, m, d):
    try:
        date(y, m, d)
        return True
    except ValueError:
        return False


def on_nfs(path="."):
    """True if `path` is on NFS, which on this cluster rejects ':' in filenames."""
    try:
        out = subprocess.run(
            ["stat", "-f", "-c", "%T", path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, timeout=5,
        ).stdout.strip()
        return out == "nfs"
    except Exception:
        return os.path.abspath(path).startswith("/cluster/VAST/")


def segments(command):
    """Split a compound shell command into argv lists."""
    try:
        toks = shlex.split(command, posix=True, comments=True)
    except ValueError:
        toks = command.split()
    segs, cur = [], []
    for t in toks:
        if t in CONTROL_TOKENS:
            if cur:
                segs.append(cur)
                cur = []
        else:
            cur.append(t)
    if cur:
        segs.append(cur)
    return segs


def strip_env_prefix(argv):
    """Drop leading VAR=value assignments and wrappers like `time`, `nohup`."""
    i = 0
    while i < len(argv):
        tok = argv[i]
        if re.match(r"^[A-Za-z_][A-Za-z0-9_]*=", tok) or tok in ("time", "nohup", "command", "exec"):
            i += 1
        elif tok == "env":
            i += 1
            while i < len(argv) and re.match(r"^[A-Za-z_][A-Za-z0-9_]*=", argv[i]):
                i += 1
        else:
            break
    return argv[i:]


def flag_value(argv, names):
    """Value of `--flag V` or `--flag=V` for any name in `names`."""
    for i, tok in enumerate(argv):
        for n in names:
            if tok == n and i + 1 < len(argv):
                return argv[i + 1]
            if tok.startswith(n + "="):
                return tok.split("=", 1)[1]
    return None


def touches_protected(argv, raw):
    for tok in argv:
        if tok.startswith("-"):
            continue
        p = os.path.abspath(os.path.expanduser(tok)) if not tok.startswith("~") else os.path.expanduser(tok)
        for pref in PROTECTED_PREFIXES:
            if p.startswith(os.path.abspath(os.path.expanduser(pref))):
                return tok
        if tok.rstrip("/") in ("results", "logs", "resources") or tok.startswith(PROTECTED_RELATIVE):
            return tok
    return None


def main():
    payload = json.load(sys.stdin)
    if payload.get("tool_name") != "Bash":
        return
    raw = (payload.get("tool_input") or {}).get("command", "") or ""
    if not raw.strip():
        return

    user = os.environ.get("USER") or os.environ.get("LOGNAME") or ""
    in_alloc = bool(os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_JOBID"))
    host = socket.gethostname()
    on_compute_host = bool(re.match(r"^[cg]\d+", host))
    on_login = not in_alloc and not on_compute_host
    try:
        cpus = int(os.environ.get("SLURM_CPUS_PER_TASK") or os.cpu_count() or 1)
    except ValueError:
        cpus = os.cpu_count() or 1

    limits = None  # lazily loaded; sinfo costs ~50ms

    for argv in segments(raw):
        argv = strip_env_prefix(argv)
        if not argv:
            continue
        prog = os.path.basename(argv[0])
        args = argv[1:]
        informational = any(a in VERSION_FLAGS for a in args)

        # --- R3: account sharing / privilege escalation ---------------------
        if prog in ("sudo", "su", "doas"):
            decide("deny", (
                f"Hellbender rule R3: `{prog}` is not permitted. Cluster policy states "
                "'Each user must use their own account to access RSS resources. Account "
                "sharing is prohibited.' There is no root access on Hellbender; if you need "
                "software installed, build a conda env under /mnt/pixstor/data/$USER or "
                "email itrss-support@umsystem.edu."
            ))
        if prog in ("ssh", "scp", "sftp", "rsync"):
            for tok in args:
                m = re.match(r"^([A-Za-z0-9._-]+)@([A-Za-z0-9._-]+)", tok)
                if not (m and user) or m.group(1) == user:
                    continue
                login, host_part = m.group(1), m.group(2)
                # Forge service accounts (git@github.com and friends) are a URL
                # convention, not a shared human login -- authentication is still
                # by the caller's own key. R3 is about people sharing UM accounts.
                if login in SERVICE_LOGINS or host_part.lower() in FORGE_HOSTS:
                    continue
                decide("deny", (
                    f"Hellbender rule R3 (no account sharing): this command connects as "
                    f"'{login}@{host_part}' but you are '{user}'. Cluster policy states "
                    "'Each user must use their own account to access RSS resources. "
                    "Account sharing is prohibited.'"
                ))

        # --- R1: no compute on the login node -------------------------------
        if on_login and prog in HEAVY_BINARIES and not informational:
            decide("deny", (
                f"Hellbender rule R1: '{prog}' is compute work and this shell is on a login "
                f"node ({host}, no SLURM allocation). ITRSS: 'Under no circumstances should "
                "your code be running on the login node.' Submit it instead:\n"
                "  sbatch --partition=general --time=<hh:mm:ss> ...\n"
                "  snakemake --profile config/slurm --use-conda <target>\n"
                "or grab an interactive allocation with `salloc` (4 h cap)."
            ))
        if on_login and prog == "python" and any("workflow/scripts/" in a for a in args):
            decide("deny", (
                "Hellbender rule R1: this runs a pipeline script on a login node. "
                "Submit it via sbatch or `snakemake --profile config/slurm`."
            ))

        # --- Snakemake must be routed through the SLURM profile -------------
        if prog == "snakemake":
            non_exec = any(a in SNAKEMAKE_NO_EXEC or a.startswith("--report=") for a in args)
            profile = flag_value(args, ["--profile"]) or ""
            executor = flag_value(args, ["--executor"]) or ""
            routed = "slurm" in profile or executor.startswith("slurm")
            if not non_exec and not routed:
                decide("deny", (
                    "Hellbender rule R1: this Snakemake invocation would execute rules "
                    "locally instead of submitting them to SLURM. Add the cluster profile:\n"
                    "  snakemake --profile config/slurm --use-conda <target>\n"
                    "Dry runs (`-n`), `--dag`, `--unlock` and `--conda-create-envs-only` are "
                    "exempt and may run as-is."
                ))
            # Conda envs must not be built on NFS: installing perl dies with
            # [Errno 22] on 'App::Cpan.3' and conda rolls the whole env back,
            # often while the shell still reports success.
            if "--use-conda" in args:
                prefix = flag_value(args, ["--conda-prefix"]) or os.environ.get(
                    "SNAKEMAKE_CONDA_PREFIX", "")
                if not prefix and on_nfs("."):
                    decide("deny", (
                        "This would build Snakemake's conda envs under .snakemake/ on the "
                        "lab VAST share, which is NFS and rejects ':' in filenames. Perl "
                        "(an orthofinder/trimal dependency) fails there with "
                        "\"[Errno 22] Invalid argument: '.../man/man3/App::Cpan.3'\" and "
                        "conda silently rolls the entire environment back. Point the envs "
                        "at pixstor (GPFS) first:\n"
                        "  export SNAKEMAKE_CONDA_PREFIX=/mnt/pixstor/data/$USER/carnivor-conda\n"
                        "or use workflow/scripts/run_snakemake.sbatch, which sets it."
                    ))

            cores = flag_value(args, ["--cores", "-c", "--local-cores"])
            if cores and cores.isdigit() and int(cores) > cpus and not routed:
                decide("deny", (
                    f"--cores {cores} exceeds this allocation's {cpus} CPUs "
                    f"(SLURM_CPUS_PER_TASK). Local execution would be throttled by the "
                    "cgroup. Route the work through --profile config/slurm instead; there "
                    "`-j` controls how many *cluster jobs* run at once, not local cores."
                ))

        # --- Thread oversubscription for locally-run tools ------------------
        if prog in THREAD_FLAGS and not informational:
            val = flag_value(args, list(THREAD_FLAGS[prog]) or ["--threads"])
            if val and val.isdigit() and int(val) > cpus:
                decide("deny", (
                    f"'{prog}' requests {val} threads but this allocation has {cpus} CPUs. "
                    f"Inside the cgroup the extra threads only add contention. Either lower "
                    f"the count to {cpus} or submit the work with sbatch/--profile "
                    "config/slurm and request more cpus-per-task there."
                ))

        # --- R2 / R7: batch submission sanity -------------------------------
        if prog in ("sbatch", "salloc", "srun"):
            if limits is None:
                limits = partition_limits()
            part = flag_value(args, ["-p", "--partition"])
            tval = flag_value(args, ["-t", "--time"])
            mins = parse_slurm_time(tval) if tval else None

            if prog == "sbatch" and part == "interactive":
                decide("deny", (
                    "Hellbender rule R2: the 'interactive' partition does not accept batch "
                    "jobs (4 h cap, interactive use only). Use --partition=general for batch "
                    "work, or `salloc --partition=interactive` for a live session."
                ))
            if part and mins and part in limits and mins > limits[part]:
                cap = limits[part]
                decide("deny", (
                    f"Hellbender rule R2: --time={tval} ({mins} min) exceeds the '{part}' "
                    f"partition cap of {cap} min ({cap/1440:.2g} days). SLURM will reject "
                    "this job. Lower the walltime, or split the work across jobs."
                ))
            if mins:
                nm = maintenance_start()
                if nm:
                    end = datetime.now() + timedelta(minutes=mins)
                    if end >= nm:
                        slack = (nm - datetime.now()).total_seconds() / 60
                        decide("ask", (
                            f"Hellbender rule R7: --time={tval} ({mins} min) would still be "
                            f"running at {nm:%Y-%m-%d %H:%M}, when the maintenance outage "
                            f"begins. Only {slack:.0f} min remain before it starts, so SLURM "
                            "will hold this job until maintenance ends rather than starting it "
                            f"now. Request at most {max(int(slack) - 5, 1)} min to run tonight, "
                            "or confirm you accept the delayed start."
                        ))

        # --- R5: irrecoverable deletion / overwrite -------------------------
        if prog in ("rm", "shred"):
            recursive = any(a.startswith("-") and ("r" in a or "R" in a) for a in args)
            target = touches_protected(args, raw)
            if target and (recursive or prog == "shred"):
                decide("ask", (
                    f"Hellbender rule R5: this would recursively delete '{target}' on cluster "
                    "storage. ITRSS: 'None of the cluster attached storage available to users "
                    "is backed up in any way by us.' There is no undo and no snapshot. "
                    "Confirm only if you are certain this data is reproducible."
                ))
        if prog == "mv":
            target = touches_protected(args, raw)
            if target and len(args) >= 2:
                pass  # moves within the tree are routine; deletion is the risk

        # --- R6: Data Class 3/4 tripwire ------------------------------------
        for tok in args:
            if tok.startswith("-"):
                continue
            if DCL_PATTERN.search(os.path.basename(tok)):
                decide("ask", (
                    f"Hellbender rule R6: '{tok}' looks like it may hold restricted data. "
                    "Cluster storage is limited to UM System Data Class 1 and 2; DCL3/DCL4 "
                    "(PHI, HIPAA, identifiable human-subject data) must not be stored here — "
                    "contact ITRSS for an approved enclave. This is a filename heuristic, so "
                    "confirm if it is a false positive."
                ))

    # Truncating redirects onto protected paths (regex on the raw string, since
    # shlex drops the redirection operator's association with its target).
    for m in re.finditer(r"(?<!\d)>(?!>)\s*([^\s;|&]+)", raw):
        tgt = m.group(1)
        if touches_protected([tgt], raw) and os.path.exists(os.path.expanduser(tgt)):
            decide("ask", (
                f"Hellbender rule R5: '>' would truncate the existing file '{tgt}' on "
                "un-backed-up cluster storage. Use '>>' to append, or confirm the overwrite."
            ))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # fail open, but loudly
        print(f"[hellbender_guard] non-fatal hook error: {exc!r}", file=sys.stderr)
    sys.exit(0)
