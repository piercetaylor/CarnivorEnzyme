# 09 — Fix silent-failure patterns in `04_run_hmmer_scan.sh`

Task: fix four issues in `workflow/scripts/braker2/04_run_hmmer_scan.sh` flagged by a prior
review pass. Real edits applied to the actual script; every fix was verified by running the
specific bash/python logic pattern in isolation with synthetic stand-in files (no real `hmmer`
binary or real proteome data available in this environment), and then by running the **actual,
edited script file** end-to-end with stubbed `hmmsearch`/`hmmfetch` binaries on `PATH`.

---

## Fix 1 (HIGH) — zero-hits guard unreachable under `set -euo pipefail`

**Root cause**: `hit_ids=$(grep -v '^#' "$hit_tbl" | ... )` — when `$hit_tbl` has zero
non-comment lines, `grep -v` exits 1, `pipefail` propagates that through the pipeline, and the
command substitution's nonzero exit kills the whole script under `set -e`, silently, before the
`[[ -z "$hit_ids" ]]` check on the next line ever runs.

**Scratch test — reproduce the bug (before fix)**:

```bash
mkdir -p /tmp/hmmer_test && cd /tmp/hmmer_test
printf '# comment line 1\n# comment line 2\n' > fake.tblout
bash -c 'set -euo pipefail; hit_ids=$(grep -v "^#" fake.tblout | awk "{print \$1}" | sort -u); echo "reached: $hit_ids"'
echo "EXIT CODE: $?"
```

**Actual output**:

```
EXIT CODE: 1
```

(No "reached:" line printed at all — confirms the script dies silently before the intended `0
hits` message, exactly as the review predicted.)

**Fix applied** (`workflow/scripts/braker2/04_run_hmmer_scan.sh`, in the per-family loop):

```bash
hit_ids=$(grep -v '^#' "$hit_tbl" | awk '{print $1}' | sort -u || true)
```

**Scratch test — confirm fix reaches the `-z` check**:

```bash
bash -c 'set -euo pipefail; hit_ids=$(grep -v "^#" fake.tblout | awk "{print \$1}" | sort -u || true); echo "reached: [$hit_ids]"'
echo "EXIT CODE: $?"
```

**Actual output**:

```
reached: []
EXIT CODE: 0
```

**Regression check** (non-empty hits still parse correctly):

```bash
printf '# comment\nSEQ1 - PF00182.20 - 1e-10 30.0 0.0 1 1 1e-09 30.0 0.0 -\nSEQ2 - PF00182.20 - 1e-08 25.0 0.0 1 1 1e-07 25.0 0.0 -\n' > fake2.tblout
bash -c 'set -euo pipefail; hit_ids=$(grep -v "^#" fake2.tblout | awk "{print \$1}" | sort -u || true); echo "reached: [$hit_ids]"'
```

Output: `reached: [SEQ1\nSEQ2]` — normal-case behavior unaffected.

---

## Fix 2 (HIGH) — all-proteomes-missing exits 0 with a misleading success summary

**Root cause**: the per-species `continue` on a missing proteome has no fallback counter; the
final `=== Summary ===` block prints unconditionally regardless of whether anything was actually
scanned.

**Scratch test — reproduce the bug (loop skeleton mirroring the real script, all species missing)**:

```bash
cat > scratch_loop.sh << 'SCRIPT'
set -euo pipefail
PROTEOME_DIR="./nonexistent_dir"
SPECIES=("sp1" "sp2" "sp3")
n_scanned=0
for species in "${SPECIES[@]}"; do
    proteome="${PROTEOME_DIR}/${species}.faa"
    if [[ ! -f "$proteome" ]]; then
        echo "[SKIP] $species — proteome not found"
        continue
    fi
    echo "=== $species ==="
    n_scanned=$((n_scanned + 1))
done
echo "n_scanned=$n_scanned"
if (( n_scanned == 0 )); then
    echo "ERROR: no proteomes found" >&2
    exit 1
fi
echo "=== Summary ==="
SCRIPT
bash scratch_loop.sh; echo "EXIT: $?"
```

**Actual output**:

```
[SKIP] sp1 — proteome not found
[SKIP] sp2 — proteome not found
[SKIP] sp3 — proteome not found
n_scanned=0
ERROR: no proteomes found
EXIT: 1
```

**Counter-placement check** (one species present, confirms `n_scanned` increments once per
successfully-scanned species, not once per attempted species regardless of skip):

```
[SKIP] sp1 — proteome not found
=== sp2 ===
[SKIP] sp3 — proteome not found
n_scanned=1
=== Summary ===
EXIT: 0
```

**Fix applied**: added `n_scanned=0` before the species loop, `n_scanned=$((n_scanned + 1))`
immediately after the `echo "=== $species ==="` line (i.e. only reached when the proteome file
exists), and the exit-1 gate immediately before the `=== Summary ===` block.

**End-to-end confirmation on the real, edited script** (all 5 species missing —
see the full-script test in Fix 3's section, second run):

```
[SKIP] darlingtonia_californica — proteome not found
[SKIP] heliamphora_ciliata — proteome not found
[SKIP] pinguicula_moranensis — proteome not found
[SKIP] drosera_spatulata — proteome not found
[SKIP] sarracenia_alata — proteome not found
ERROR: no proteomes found under resources/accessions/tarnita2023 — run 00_download_genomes.sh through 03_extract_proteins.sh first
EXIT CODE: 1
```

---

## Fix 3 (MEDIUM) — unresolvable hmmsearch hit IDs silently vanish, message blames length filter

**Root cause**: the embedded Python does
`hits = [records[hid] for hid in hit_ids if hid in records]` (silently drops any `.tblout` hit ID
absent from the parsed proteome FASTA — e.g. stale `.tblout` from a different proteome file),
then the printed summary only ever reported `len(hit_ids)` and `len(kept)`, so IDs dropped at the
existence-filter step were invisible and any shortfall looked like it came from the length filter.

**Scratch test — reproduce the bug (old logic, `fake_proteome.faa` has SEQ1 (long) + SEQ2
(short); `hit_ids` includes a third ID, SEQ3, that isn't in the FASTA at all)**:

```python
from Bio import SeqIO
proteome = "fake_proteome.faa"
hit_ids = set("SEQ1 SEQ2 SEQ3".strip().split())
min_len = 200
records = {r.id: r for r in SeqIO.parse(proteome, "fasta")}
hits = [records[hid] for hid in hit_ids if hid in records]
kept = [r for r in hits if len(r.seq) >= min_len]
print(f"  family: {len(hit_ids)} hmmer hits -> {len(kept)} kept (>={min_len} aa) -> out_old.faa")
```

**Actual output**:

```
  family: 3 hmmer hits -> 1 kept (>=200 aa) -> out_old.faa
```

This message implies 2 of 3 hits were removed by the length cutoff. In reality only SEQ2 (short)
was length-filtered; SEQ3 was silently dropped because it doesn't exist in the FASTA (simulating a
stale `.tblout`) — exactly the misattribution the review flagged.

**Fix applied**:

```python
if len(hits) < len(hit_ids):
    print(f"  WARNING: {len(hit_ids) - len(hits)} hmmer hit ID(s) not found in proteome FASTA — check .tblout matches the proteome file")
print(f"  ${family}: {len(hit_ids)} hmmer hits → {len(hits)} resolved → {len(kept)} kept (≥{min_len} aa) → {out_fa}")
```

**Scratch test — confirm fix (new logic, same synthetic inputs)**:

```
  WARNING: 1 hmmer hit ID(s) not found in proteome FASTA — check .tblout matches the proteome file
  test_family: 3 hmmer hits -> 2 resolved -> 1 kept (>=200 aa) -> out_new.faa
```

Now correctly shows 3 hits → 2 resolved (1 unresolvable, warned) → 1 kept (by length).

**End-to-end confirmation** — ran the actual, edited `04_run_hmmer_scan.sh` with stubbed
`hmmfetch`/`hmmsearch` on `PATH` (`/tmp/hmmer_full/bin`), a synthetic
`resources/accessions/tarnita2023/darlingtonia_californica.faa` (one 401-aa "SEQ_LONG", one 51-aa
"SEQ_SHORT"), and a stub `hmmsearch` that emits different canned `.tblout` content per Pfam
accession requested (zero hits for PF00182 → tests Fix 1; a hit plus one stale/nonexistent ID for
PF00149 → tests Fix 3; a normal single resolvable hit for PF00445; a hit that resolves but fails
the length filter for PF00026):

```
=== darlingtonia_californica ===
hmmsearch stub ran for resources/hmm/PF00026.hmm vs resources/accessions/tarnita2023/darlingtonia_californica.faa
  aspartic_proteases_a1b_homology: 1 hmmer hits → 1 resolved → 0 kept (≥280 aa) → resources/accessions/tarnita2023_by_family/darlingtonia_californica_aspartic_proteases_a1b_homology.faa
hmmsearch stub ran for resources/hmm/PF00445.hmm vs resources/accessions/tarnita2023/darlingtonia_californica.faa
  rnase_t2: 1 hmmer hits → 1 resolved → 1 kept (≥180 aa) → resources/accessions/tarnita2023_by_family/darlingtonia_californica_rnase_t2.faa
hmmsearch stub ran for resources/hmm/PF00182.hmm vs resources/accessions/tarnita2023/darlingtonia_californica.faa
  chitinases_gh19_class_iv: 0 hits
hmmsearch stub ran for resources/hmm/PF00149.hmm vs resources/accessions/tarnita2023/darlingtonia_californica.faa
  WARNING: 1 hmmer hit ID(s) not found in proteome FASTA — check .tblout matches the proteome file
  purple_acid_phosphatase: 2 hmmer hits → 1 resolved → 1 kept (≥350 aa) → resources/accessions/tarnita2023_by_family/darlingtonia_californica_purple_acid_phosphatase.faa
[SKIP] heliamphora_ciliata — proteome not found
[SKIP] pinguicula_moranensis — proteome not found
[SKIP] drosera_spatulata — proteome not found
[SKIP] sarracenia_alata — proteome not found

=== Summary ===
Family hit FASTAs written to resources/accessions/tarnita2023_by_family
...
EXIT CODE: 0
```

This single run exercises Fix 1 (`chitinases_gh19_class_iv: 0 hits` — no silent death), Fix 2
(one species present → `n_scanned=1` → normal exit 0 summary), and Fix 3 (`purple_acid_phosphatase`
WARNING + correct 3-count message) simultaneously, against the real edited file.

Removing the proteome entirely and re-running (all 5 species skipped) reproduces the Fix 2 exit-1
path shown above.

**Aside (not part of the assigned fixes, noted for completeness)**: the first attempt at this
end-to-end run hit an unrelated pre-existing environment issue — Python on this Windows box
defaults `stdout` to `cp1252`, which can't encode the `→`/`≥` characters already used in the
script's f-strings, raising `UnicodeEncodeError`. This bug pre-dates all four fixes (the arrow
character was already present in the original `print` statement) and is orthogonal to the task,
so it was worked around for testing purposes only (`PYTHONIOENCODING=utf-8`) and left unfixed in
the script, since it wasn't in scope. Flagging here in case the PI wants a follow-up task for it —
likely relevant for any Windows-side manual runs of this script (its intended runtime is Linux
HPC, per the header comment `Run from: /home/pmt5gt/data/CarnivorEnzyme`).

---

## Fix 4 (LOW) — stale example path in the final "Next steps" echo block

**Root cause**: line 137 (pre-fix) still referenced the pre-split family name
(`darlingtonia_californica_chitinases_gh19.faa`) in the printed example, even though the
associative-array keys were already renamed to `chitinases_gh19_class_iv` in an earlier pass (per
`audit/CHANGELOG.md` T3 entry).

**Verification**: grepped the file for the stale bare pattern after editing:

```bash
grep -nE 'chitinases_gh19[^_]' "workflow/scripts/braker2/04_run_hmmer_scan.sh"
```

**Actual output**: `No matches found` (confirms zero remaining stale references anywhere in the
file).

**Fix applied**:

```
echo "             - local:resources/accessions/tarnita2023_by_family/darlingtonia_californica_chitinases_gh19_class_iv.faa"
```

**Duplicate-comment check**: the manual Class I vs IV triage note already exists in the
`FAMILY_HMMS` comment block (lines 35-38: "a fresh hmmer hit still needs manual Class I vs IV
domain-architecture triage ... before being filed under chitinases_gh19_class_iv vs
chitinases_gh19_class_i — this script only finds candidates, it doesn't classify them"). No
duplicate note was added elsewhere.

---

## Syntax check

```bash
bash -n "workflow/scripts/braker2/04_run_hmmer_scan.sh" && echo "SYNTAX OK"
```

Output: `SYNTAX OK`
