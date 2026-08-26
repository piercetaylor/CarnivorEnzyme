# 13 — Distinguish real I/O errors from legitimate zero-hits in `04_run_hmmer_scan.sh`

## Background

`audit/09_hmmer_scan_script_fixes.md` added `|| true` to the hit-ID extraction pipeline:

```bash
hit_ids=$(grep -v '^#' "$hit_tbl" | awk '{print $1}' | sort -u || true)
```

This was needed because under `set -euo pipefail`, `grep -v '^#'` exits `1` when a `.tblout`
file has zero non-comment lines — a legitimate "hmmsearch found nothing for this family in this
species" outcome — and an unguarded pipeline failure there would kill the whole script before the
`-z` check for "0 hits" ever ran.

Follow-up review found the fix over-broadened: `|| true` on the whole pipeline swallows **any**
nonzero exit from **any** stage of `grep | awk | sort` under `pipefail`, including:

- `grep` exit code `2` (file missing, unreadable, or a real I/O error) — semantically different
  from exit code `1` ("no match"), but `|| true` can't distinguish them.
- A failure in `awk` or `sort` themselves.

All of these were previously reported identically as `"$family: 0 hits"` and the loop continued,
meaning a real infrastructure/permissions/disk problem was indistinguishable from a legitimate
empty hmmer result.

## Fix

Replaced the single guarded pipeline with an isolated `grep` call whose exit code is inspected
before deciding whether to continue or abort:

```bash
set +e
raw_hits=$(grep -v '^#' "$hit_tbl")
grep_rc=$?
set -e
if (( grep_rc > 1 )); then
    echo "ERROR: failed to read $hit_tbl (grep exit code $grep_rc)" >&2
    exit 1
fi
hit_ids=$(printf '%s\n' "$raw_hits" | awk '{print $1}' | sort -u)
```

- `grep_rc == 0` (matches found) or `grep_rc == 1` (no matches — legitimate empty result) both
  fall through to the `awk | sort` step as before.
- `grep_rc >= 2` (real error: missing file, unreadable file, etc.) now aborts the script loudly
  with a descriptive message instead of being silently reported as "0 hits".
- `awk`/`sort` are no longer wrapped in a blanket `|| true` — if the file was readable, these are
  expected to succeed on well-formed input, and failures there will now propagate under
  `set -euo pipefail` instead of being masked.
- `printf '%s\n' "$raw_hits" | awk ... | sort -u` on an empty `$raw_hits` still correctly resolves
  to an empty `hit_ids` (verified in Test 1 below), so the existing `[[ -z "$hit_ids" ]]` "0 hits,
  continue" branch is preserved for the legitimate empty case.

## Verification

All three scratch tests were run from a temp scaffold script (`/tmp/hmmer_scratch/test_snippet.sh`)
that reproduces the exact fixed logic block in isolation, run against synthetic `.tblout`-like
files — no real repo data was touched.

Scratch harness:

```bash
mkdir -p /tmp/hmmer_scratch && cd /tmp/hmmer_scratch

cat > test_snippet.sh <<'SCRIPT'
#!/usr/bin/env bash
set -euo pipefail

hit_tbl="$1"

set +e
raw_hits=$(grep -v '^#' "$hit_tbl")
grep_rc=$?
set -e
echo "grep_rc=$grep_rc"
if (( grep_rc > 1 )); then
    echo "ERROR: failed to read $hit_tbl (grep exit code $grep_rc)" >&2
    exit 1
fi
hit_ids=$(printf '%s\n' "$raw_hits" | awk '{print $1}' | sort -u)
echo "hit_ids=[${hit_ids}]"
if [[ -z "$hit_ids" ]]; then
    echo "RESULT: 0 hits (tolerated, continuing)"
else
    echo "RESULT: hits found"
fi
SCRIPT
chmod +x test_snippet.sh
```

### Test 1 — comment-only `.tblout` (legitimate zero hits, `grep` exit 1)

```bash
cat > comment_only.tbl <<'EOF'
# This is a comment
# target name        accession  query name
#------------------- ---------- -----------
EOF
./test_snippet.sh comment_only.tbl
echo "script exit code: $?"
```

Output:

```
grep_rc=1
hit_ids=[]
RESULT: 0 hits (tolerated, continuing)
script exit code: 0
```

Confirms: `grep_rc == 1` is correctly tolerated, `hit_ids` resolves to an empty string (not a
single blank-line artifact — `printf '%s\n' ""` piped through `awk '{print $1}' | sort -u`
collapses to nothing once the trailing newline is stripped by command substitution), and the
script exits `0` instead of aborting.

### Test 2 — nonexistent file (real I/O error, `grep` exit 2)

```bash
set +e
./test_snippet.sh /tmp/hmmer_scratch/does_not_exist.tbl
echo "script exit code: $?"
set -e
```

Output:

```
grep: /tmp/hmmer_scratch/does_not_exist.tbl: No such file or directory
grep_rc=2
ERROR: failed to read /tmp/hmmer_scratch/does_not_exist.tbl (grep exit code 2)
script exit code: 1
```

Confirms: `grep_rc == 2` is now caught and reported with the new `ERROR:` message, and the script
exits `1` instead of silently reporting "0 hits" as it would have under the old blanket `|| true`.

### Test 3 — normal `.tblout` with real hit lines

```bash
cat > real_hits.tbl <<'EOF'
# target name        accession  query name
#------------------- ---------- -----------
protein_A            -          PF00182   -    1.2e-40  130.5   0.2  1  1  1.5e-42  9e-40  135.0   0.1     1   200     5   210     3   215 0.97 -
protein_B            -          PF00182   -    3.4e-50  160.1   0.0  1  1  4.0e-51  2e-49  159.0   0.0     1   200     1   198     1   200 0.99 -
protein_A            -          PF00182   -    5.0e-38  120.0   0.0  1  1  6.0e-38  1e-37  118.0   0.0     1   190     8   200     5   205 0.95 -
EOF
./test_snippet.sh real_hits.tbl
echo "script exit code: $?"
```

Output:

```
grep_rc=0
hit_ids=[protein_A
protein_B]
RESULT: hits found
script exit code: 0
```

Confirms: the normal case (including a duplicate ID across two hit lines) still correctly
resolves to a deduplicated, sorted `hit_ids` list (`protein_A`, `protein_B`), matching prior
behavior.

## Conclusion

The three scratch tests confirm the fix correctly distinguishes:

1. Legitimate zero-hit `.tblout` files (`grep` exit `1`) — tolerated, loop continues.
2. Real I/O errors such as a missing/unreadable file (`grep` exit `2`) — now caught and aborted
   with a descriptive `ERROR:` message instead of being misreported as "0 hits".
3. The normal case with real hit lines — unaffected, same output as before.

Applied to `workflow/scripts/braker2/04_run_hmmer_scan.sh` (hit-ID extraction block, inside the
per-family loop).
