# 12 — Make silently-empty species/families visible in `fetch_sequences.py`

Follow-up to `audit/08_fetch_sequences_hardening.md` (T8), which fixed five silent-failure
patterns in `workflow/scripts/fetch_sequences.py` but left two more: a species whose every fetched
record gets filtered out by the length cutoff, and a family whose every accession is `TODO`,
both wrote no output and left no trace in the exit code. Both are fixed below, each verified by
actually running the real CLI (never just read-and-reason) — once against the real
`config/enzyme_families.yaml` with live NCBI/UniProt fetches, and once against a disposable
scratch YAML built to force both conditions deterministically.

---

## What the real config actually contains

Before writing any code, checked whether any *currently-configured* family (`tier1` +
`methods_benchmark` — the two sections `fetch_sequences.py` actually iterates; `tier2` is defined
in the YAML but never fetched) has zero non-TODO accessions in total, which would make Fix 2 an
urgent, live-firing bug rather than a prophylactic guard.

Counted TODO vs. non-TODO accession lines directly against the real file:

```
$ python3 -c "
import yaml
config = yaml.safe_load(open('config/enzyme_families.yaml', encoding='utf-8'))
for section in ('tier1', 'methods_benchmark'):
    for fam, data in config.get(section, {}).items():
        total = todo = 0
        for accs in data.get('accessions', {}).values():
            for a in accs:
                total += 1
                if not a or str(a).strip().upper().startswith('TODO'):
                    todo += 1
        print(f'{section}.{fam}: {total} accessions, {todo} TODO, {total-todo} usable')
"
tier1.chitinases_gh19_class_iv: 27 accessions, 10 TODO, 17 usable
tier1.purple_acid_phosphatase: 21 accessions, 5 TODO, 16 usable
tier1.rnase_t2: 22 accessions, 6 TODO, 16 usable
tier1.aspartic_proteases_a1b_homology: 21 accessions, 4 TODO, 17 usable
methods_benchmark.chitinases_gh19_class_i: 3 accessions, 0 TODO, 3 usable
methods_benchmark.nepenthesins: 14 accessions, 0 TODO, 14 usable
methods_benchmark.neprosins: 8 accessions, 3 TODO, 5 usable
```

**Finding: none of the 7 currently-fetched families has zero non-TODO accessions.** Every family
has at least 3 usable accessions. Fix 2 is therefore **prophylactic** today (guards against a
future edit that TODOs-out an entire family, or a `tier2` family being wired into the fetch list
without checking), not an active bug — confirmed by direct count, not by inspection alone.

---

## Fix 1: species that fetch successfully but get entirely length-filtered wrote nothing, invisibly

**Before** (around line 304): `if not kept: continue` after the per-species length filter. If
every record for a species survived Pass 1 (fetched fine) but got excluded in Pass 2 (all below
`min_length_fraction * median`), the species silently vanished — no file, and while `n_skipped`
incremented per filtered *record*, nothing distinguished "this species has zero surviving records"
from "this species had one long record and one short record filtered, one file still written."
`n_skipped` also never affected the exit code (only `n_failed` did), so this was invisible in both
the log-scanning sense and the automation/exit-code sense.

**Judgment call — should this also force `n_failed`/`sys.exit(1)`?** No. Several accessions
already in `config/enzyme_families.yaml` are documented, expected partial/truncated gene models
that this exact filter is designed to catch — e.g. `aspartic_proteases_a1b_homology`'s
`Dionaea_muscipula` entry: *"NOTE: 300 aa is at the lower edge of expected range [280, 650]. Very
likely a partial/truncated gene model. Flag for length filter in pipeline."* Forcing a hard
pipeline failure every time a documented partial isoform gets filtered would make the script
un-runnable against the current config (see live evidence below — this exact case fires on a real
run) and would punish the config authors for doing exactly what the length filter is for. Instead,
this is surfaced with a `WARNING`-level line at the point of occurrence plus a distinct,
impossible-to-miss aggregate counter in the final run summary, satisfying "prefer visibility over
silence" without turning an expected, already-documented outcome into a hard failure.

**After** (`workflow/scripts/fetch_sequences.py`, in the Pass-2 species loop, ~line 298-328):

```python
if not kept:
    n_species_all_filtered += 1
    logger.warning(
        "ZERO SURVIVING %s / %s: %d/%d fetched record(s) below length cutoff "
        "(%.0f aa) — no FASTA written for this species",
        family_key, species, len(records), len(records), cutoff,
    )
    continue
```

plus a new counter `n_species_all_filtered` (initialized alongside `n_fetched`/`n_skipped`/
`n_failed`) that is reported in the final summary (see below).

---

## Fix 2: a family with every accession still TODO logged one WARNING and moved on invisibly

**Before** (~line 280-282): `if not all_lengths: logger.warning(...); continue`. A single WARNING
line, no counter, no summary-level visibility, `n_failed` untouched — indistinguishable at exit
time from a fully successful run.

**After**: the per-family message is elevated to `ERROR` (a family producing literally zero
output across every species is a more severe condition than an individual filtered record), and
the family key is recorded in a new `families_zero_output: list[str]` collected across the whole
run:

```python
if not all_lengths:
    logger.error(
        "ZERO OUTPUT family %s — no non-TODO accession resolved to a sequence "
        "(check YAML TODOs); no FASTA written for any species in this family",
        family_key,
    )
    families_zero_output.append(family_key)
    continue
```

**Extended beyond the original ask:** while implementing this, noticed the same "zero output"
condition can also arise *after* Pass 1 succeeds — if a family fetches records for every species
but every single species then gets entirely length-filtered in Pass 2 (Fix 1's condition, but for
*every* species in the family simultaneously). That would leave an empty `results/sequences/
{family}/` directory with only the original single WARNING line (now per-species WARNINGs) and no
family-level signal at all — the exact same invisibility Fix 2 targets, just reached via a
different path. Added a `family_wrote_any` flag, set `True` on every successful `SeqIO.write`,
checked after the per-species loop:

```python
if not family_wrote_any:
    logger.error(
        "ZERO OUTPUT family %s — every species had all records filtered by length; "
        "no FASTA written for any species in this family",
        family_key,
    )
    families_zero_output.append(family_key)
```

Neither branch forces `sys.exit(1)` by itself (same reasoning as Fix 1 — a family being reduced to
zero output by the length filter is a length-QC outcome, not necessarily a fetch failure), but the
new `families_zero_output` list is unconditionally reported at `ERROR` level in the final summary,
so it can never be missed by someone reading a truncated log (grep for `ERROR` or `ZERO OUTPUT`
finds it immediately; the aggregate line names every affected family key by name).

---

## Final summary changes

`Done. fetched=%d skipped/filtered=%d failed=%d` became:

```python
logger.info(
    "Done. fetched=%d  skipped/filtered=%d  failed=%d  species_all_filtered=%d  "
    "families_zero_output=%d",
    n_fetched, n_skipped, n_failed, n_species_all_filtered, len(families_zero_output),
)
if n_species_all_filtered:
    logger.warning(
        "%d species had ZERO surviving sequences after length filtering (fetched "
        "successfully but every record was below the length cutoff) — see 'ZERO "
        "SURVIVING' lines above for which species/family",
        n_species_all_filtered,
    )
if families_zero_output:
    logger.error(
        "%d famil(y/ies) produced ZERO output sequences: %s",
        len(families_zero_output), families_zero_output,
    )
if n_failed > 0:
    ...  # unchanged — still the only thing that forces sys.exit(1)
```

Both new counters are now present in every run's terminal summary line regardless of `--verbose`,
and `families_zero_output` names the specific family keys so a human doesn't have to scroll back
through the log to find which ones.

---

## Verification (a): real run against the real config

Ran the actual CLI against the real `config/enzyme_families.yaml`, all 7 currently-fetched
families, real NCBI Entrez + UniProt REST fetches, credentials from `.env`, scratch output
directory (`/tmp/scratch_fetch_test12`, deleted after):

```
$ python3 workflow/scripts/fetch_sequences.py \
    -f config/enzyme_families.yaml \
    -o /tmp/scratch_fetch_test12 \
    -e "$NCBI_EMAIL" -k "$NCBI_API_KEY"
```

**Real output (relevant lines):**

```
WARNING __main__: FILTERED   Sarracenia_purpurea|KAL6998072.1 (428 aa) < 80% of median (539 aa)
WARNING __main__: ZERO SURVIVING purple_acid_phosphatase / Sarracenia_purpurea: 1/1 fetched record(s) below length cutoff (431 aa) — no FASTA written for this species
...
WARNING __main__: FILTERED   Dionaea_muscipula|A0A0E3GLN3 (352 aa) < 80% of median (456 aa)
WARNING __main__: FILTERED   Dionaea_muscipula|A0A0E3M338 (300 aa) < 80% of median (456 aa)
WARNING __main__: ZERO SURVIVING aspartic_proteases_a1b_homology / Dionaea_muscipula: 2/2 fetched record(s) below length cutoff (365 aa) — no FASTA written for this species
...
INFO __main__: Done. fetched=92  skipped/filtered=30  failed=0  species_all_filtered=2  families_zero_output=0
WARNING __main__: 2 species had ZERO surviving sequences after length filtering (fetched successfully but every record was below the length cutoff) — see 'ZERO SURVIVING' lines above for which species/family
EXIT CODE: 0
```

**This confirms two real, previously-invisible occurrences of the Fix-1 condition** in the current
config (not just a hypothetical the fix guards against):

1. `purple_acid_phosphatase / Sarracenia_purpurea` — `KAL6998072.1` (428 aa) falls just under
   that family's 431 aa cutoff. This is the exact accession/family the T5 checker
   (`audit/05_live_refetch_verification.md`) already found dropped by the filter during count
   reconciliation, but that prior pass only surfaced it via manual arithmetic after the fact — with
   this fix it is now visible directly in the run's own log output.
2. `aspartic_proteases_a1b_homology / Dionaea_muscipula` — **both** of that species' two
   accessions (Dionain-1, 352 aa; Dionain-3, 300 aa) fall under the family's 365 aa cutoff,
   entirely removing `Dionaea_muscipula` from this family's output. This is exactly the outcome
   the config's own inline comment on `A0A0E3M338` anticipated ("Very likely a partial/truncated
   gene model. Flag for length filter in pipeline") — now it actually is flagged, at runtime, not
   just in a code comment.

`families_zero_output=0` confirms the earlier count (no currently-configured family has zero
non-TODO accessions) holds true end-to-end through a live fetch, not just at the YAML-counting
level.

`fetched=92 skipped/filtered=30 failed=0` matches the T5 baseline run exactly (same live data,
same filter arithmetic) — this change is purely additive visibility, not a behavior change to
what gets fetched, filtered, or written.

## Verification (b): synthetic scratch YAML forcing both conditions

Built a disposable YAML (`scratch_test12_visibility.yaml`, deleted after the test — never touched
`config/enzyme_families.yaml`) with one family designed to hit Fix 1 (one long real sequence, one
real-but-tiny fragment) and one family designed to hit Fix 2 (both species entirely `TODO`):

```yaml
tier1:
  test_family_short_filter:
    expected_length_aa: [100, 500]
    accessions:
      Good_species:
        - BAD07474.1    # real, 437 aa
      Short_species:
        - CAA28969.1    # real, 27 aa fragment — should get length-filtered out entirely
  test_family_all_todo:
    expected_length_aa: [100, 500]
    accessions:
      Species_one:
        - TODO
      Species_two:
        - TODO
```

```
$ python3 workflow/scripts/fetch_sequences.py \
    -f scratch_test12_visibility.yaml \
    -o /tmp/scratch_fetch_test12_synth \
    -e "$NCBI_EMAIL" -k "$NCBI_API_KEY" --verbose
```

**Real output (full, both families):**

```
INFO __main__: ━━━ Family: test_family_short_filter (expected 100–500 aa) ━━━
DEBUG __main__: → NCBI: BAD07474.1
INFO __main__:   FETCHED  Good_species | BAD07474.1 (437 aa)
DEBUG __main__: → NCBI: CAA28969.1
WARNING __main__: VERY SHORT CAA28969.1: 27 aa (expected 100–500) — included but flagged
INFO __main__:   FETCHED  Short_species | CAA28969.1 (27 aa)
INFO __main__: Family test_family_short_filter: 2 sequences, median 232 aa, length cutoff 186 aa
INFO __main__:   WROTE    Good_species → ...\test_family_short_filter\Good_species.fa (1 record(s))
WARNING __main__: FILTERED   Short_species|CAA28969.1 (27 aa) < 80% of median (232 aa)
WARNING __main__: ZERO SURVIVING test_family_short_filter / Short_species: 1/1 fetched record(s) below length cutoff (186 aa) — no FASTA written for this species
INFO __main__: ━━━ Family: test_family_all_todo (expected 100–500 aa) ━━━
DEBUG __main__: SKIP TODO  test_family_all_todo / Species_one
DEBUG __main__: SKIP TODO  test_family_all_todo / Species_two
ERROR __main__: ZERO OUTPUT family test_family_all_todo — no non-TODO accession resolved to a sequence (check YAML TODOs); no FASTA written for any species in this family
INFO __main__: Done. fetched=2  skipped/filtered=3  failed=0  species_all_filtered=1  families_zero_output=1
WARNING __main__: 1 species had ZERO surviving sequences after length filtering (fetched successfully but every record was below the length cutoff) — see 'ZERO SURVIVING' lines above for which species/family
ERROR __main__: 1 famil(y/ies) produced ZERO output sequences: ['test_family_all_todo']
EXIT CODE: 0
```

Both new signals fire exactly as designed: `species_all_filtered=1` / the `ZERO SURVIVING` line
for `Short_species` (Fix 1), and `families_zero_output=1` / the `ZERO OUTPUT` line naming
`test_family_all_todo` (Fix 2). Exit code stayed `0` since nothing actually *failed* to fetch
(`failed=0`) — consistent with the design decision above that these are visibility signals, not
hard failures, given both conditions can be legitimate/expected outcomes.

---

## Cleanup

`/tmp/scratch_fetch_test12/`, `/tmp/scratch_fetch_test12_synth/`, and
`scratch_test12_visibility.yaml` were all deleted after verification. `results/sequences/` (the
real pipeline output) was never touched by any command in this task — every run above used a
scratch `-o` output directory.

---

## Summary of execution evidence

| Fix | Real command run | Real result |
|-----|------------------|-------------|
| Config TODO/usable count | Direct PyYAML introspection of all 7 fetched families | 0 of 7 families have zero non-TODO accessions — Fix 2 is prophylactic today |
| Fix 1 (species all-filtered) | Real CLI, real config, live NCBI/UniProt fetch, scratch output dir | 2 real occurrences found and surfaced (`purple_acid_phosphatase/Sarracenia_purpurea`, `aspartic_proteases_a1b_homology/Dionaea_muscipula`); `species_all_filtered=2` in summary; exit 0 |
| Fix 2 (family zero output) | Same real run as above | `families_zero_output=0` — confirms no live family currently hits this; extended fix (empty-after-filter case) did not fire either |
| Both fixes, synthetic | Real CLI against disposable scratch YAML forcing both conditions | `species_all_filtered=1`, `families_zero_output=1`, both named explicitly in ERROR/WARNING summary lines; exit 0 |
