#!/usr/bin/env python3
"""Behaviour tests for .claude/hooks/hellbender_guard.py.

Run from the repo root:  python3 .claude/hooks/test_hellbender_guard.py

Each case pipes a synthetic PreToolUse payload at the hook and asserts the
decision. "none" means the hook stayed silent and the harness applies its normal
permission logic.

Two cases are time-sensitive and assume a maintenance reservation exists:
"job running into maintenance" expects `ask`, and "short sbatch finishing
tonight" expects `none`. Immediately after a maintenance window closes, the next
reservation may not be posted yet -- if those two flip, check
`scontrol show reservation` before assuming the hook regressed.
"""

import json, os, subprocess, sys
HOOK = ".claude/hooks/hellbender_guard.py"
CASES = [
    ("deny", "sbatch 3 days on general (R2)",        "sbatch --partition=general --time=3-00:00:00 x.sh"),
    ("deny", "sbatch on interactive partition (R2)", "sbatch --partition=interactive --time=01:00:00 x.sh"),
    ("deny", "gpu partition over 2-day cap (R2)",    "sbatch --partition=gpu --time=5-00:00:00 x.sh"),
    ("deny", "unrouted snakemake (R1)",              "snakemake phase2 --cores 64"),
    ("deny", "orthofinder 128 threads",              "orthofinder -t 128 -f proteomes/"),
    ("deny", "iqtree -T 64 in a 4-cpu alloc",        "iqtree3 -s aln.fa -T 64"),
    ("deny", "ssh as another user (R3)",             "ssh someoneelse@hellbender.rnet.missouri.edu"),
    ("deny", "sudo (R3)",                            "sudo dnf install iqtree"),
    ("ask",  "rm -rf results (R5)",                  "rm -rf results/convergence"),
    ("ask",  "rm -rf on VAST abs path (R5)",         "rm -rf /cluster/VAST/mendozacozatld-lab/PierceTaylor/x"),
    ("ask",  "DCL tripwire (R6)",                    "cp secret_phi_data.csv results/"),
    ("ask",  "job running into maintenance (R7)",    "sbatch --partition=general --time=1-00:00:00 x.sh"),
    ("none", "dry run",                              "snakemake -n phase2"),
    ("none", "squeue",                               "squeue -u $USER"),
    ("none", "sinfo",                                "sinfo -s"),
    ("none", "short sbatch finishing tonight",       "sbatch --partition=general --time=03:00:00 x.sh"),
    ("none", "the real T0.1 command",                "snakemake --profile config/slurm --use-conda -j 20 results/convergence/chitinases_gh19_class_iv.convergent_sites.tsv", True),
    ("deny", "T0.1 command without the prefix set",  "snakemake --profile config/slurm --use-conda -j 20 results/convergence/chitinases_gh19_class_iv.convergent_sites.tsv", False),
    ("none", "the --touch step",                     "snakemake --profile config/slurm --use-conda --touch results/convergence/x.tsv", True),
    # --use-conda without a conda-prefix on the NFS share is the 15-minute trap
    # (perl's 'App::Cpan.3' hits [Errno 22] and conda rolls the env back).
    ("deny", "use-conda on NFS w/o prefix",          "snakemake --conda-create-envs-only --use-conda results/x.tsv", False),
    ("none", "use-conda with explicit prefix",       "snakemake --conda-create-envs-only --use-conda --conda-prefix /mnt/pixstor/data/x results/x.tsv", False),
    ("none", "iqtree3 --version",                    "iqtree3 --version"),
    ("none", "phylogeny path is not a PHI match",    "ls results/phylogenies/x.asr.state"),
    ("none", "git@github.com is a service login",     "ssh -T git@github.com"),
    ("none", "git push over an SSH remote",           "git push git@github.com:piercetaylor/CarnivorEnzyme.git main"),
    ("deny", "a real shared UM account",              "ssh otherstudent@hellbender-login.rnet.missouri.edu"),
    ("none", "rsync as myself",                      "rsync -av pmt5gt@hellbender-dtn-p1.rnet.missouri.edu:/x ."),
    ("none", "compound cmd, all clean",              "cd results && ls -la && squeue -u $USER"),
    ("deny", "compound cmd hiding a violation",      "cd results && snakemake phase2 --cores 64"),
]
PREFIX = "/mnt/pixstor/data/%s/carnivor-conda" % (os.environ.get("USER", "user"),)

fails = 0
for case in CASES:
    expected, label, cmd = case[0], case[1], case[2]
    needs_prefix = case[3] if len(case) > 3 else False
    env = dict(os.environ)
    env.pop("SNAKEMAKE_CONDA_PREFIX", None)
    if needs_prefix:
        env["SNAKEMAKE_CONDA_PREFIX"] = PREFIX
    payload = json.dumps({"tool_name": "Bash", "tool_input": {"command": cmd}})
    proc = subprocess.run([sys.executable, HOOK], input=payload, env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True)
    out = proc.stdout.strip()
    got = json.loads(out)["hookSpecificOutput"]["permissionDecision"] if out else "none"
    ok = got == expected
    fails += not ok
    print("  %s  exp=%-5s got=%-5s  %s" % ("PASS" if ok else "FAIL", expected, got, label))
    if proc.stderr.strip():
        print("        stderr: " + proc.stderr.strip().splitlines()[0])
print("\n%d/%d passed" % (len(CASES) - fails, len(CASES)))
sys.exit(1 if fails else 0)
