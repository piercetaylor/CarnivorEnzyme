# Audit 04: Split GH19 chitinases into Class IV (convergence target) vs. Class I (structural only)

**Date:** 2026-08-22
**Task:** T4 — verify the Class I vs Class IV domain-architecture claim for the `chitinases_gh19`
tier1 entry in `config/enzyme_families.yaml`, and, gated on that verification, split the entry so
Class I sequences are no longer mixed into the Class IV cross-lineage convergence alignment.

---

## Task 1: Verify the Class I vs Class IV domain-architecture claim

### Source 1 — Yoneda et al. 2025, *FEBS Open Bio* 15:1930–1944 (PDB 9JTR)

Fetched directly: https://febs.onlinelibrary.wiley.com/doi/full/10.1002/2211-5463.70110

> "Plant chitinases are further divided into five major classes based on the similarities in
> their primary amino acid sequences. Sequence analysis of *D. adelae* chitinase has revealed
> that it belongs to class I of the GH19 family. A distinguishing feature of class I chitinases
> is the presence of a chitin-binding domain, also known as the hevein domain, at their
> N-terminal region."

> "Proteomic analysis has revealed the presence of a class I chitinase belonging to the
> glycosyl hydrolase 19 family (GH19) in the digestive fluid of the carnivorous plant *D. adelae*."

The paper identifies the protein as GenBank **BAW35424.1** and confirms "class I" throughout, with
the crystal structure (catalytic domain) deposited as PDB 9JTR (1.73 Å / 1.57 Å datasets).

This directly confirms: (a) Class I is genuinely defined by an intact N-terminal chitin-binding
(hevein/CBM18) domain, and (b) D. adelae's BAW35424.1 is Class I, not Class IV — exactly matching
what the config file's own prior comments already claimed (`class_i=true`).

### Source 2 — GH19 Engineering Database (Karlsson & Stenlid framework; PMC8547705 / PLOS ONE)

WebSearch of the domain-architecture literature returned a direct definitional quote (paraphrased
by the search tool from the primary source, consistent with the paper's known content):

> "Class I chitinase consists of a CBM18 and 'loopful' GH19 catalytic domain; class II chitinase
> consists of a loopful GH19 catalytic domain only; class IV chitinase consists of a CBM18 with a
> deletion and 'loopless' GH19 catalytic domain; class II-L chitinase consists of 'loopless' GH19
> catalytic domain only."

This confirms the config's "truncated/absent CBD" framing for Class IV: Class IV retains a
CBM18 *with a deletion* (truncated), distinct from Class I's intact CBM18 — and confirms Class I
and Class IV differ in both the chitin-binding domain and the catalytic-domain loop structure
("loopful" vs. "loopless"), i.e. two independent architectural differences, not one cosmetic label.

### Source 3 — Drosera binata DbChitI-3 (Planta 2025 / PMC11725546)

WebSearch of "Complex transcription regulation of acidic chitinase suggests fine-tuning of
digestive processes in *Drosera binata*" (Planta 2025, PMC11725546) returned:

> "DbChitI-3 is a novel chitinase gene isolated from the carnivorous plant species *Drosera
> binata* with strong homology to other *Drosera* species' extracellular class I chitinases...
> The protein is identified with accession number QHN63861.1."

This independently confirms D. binata's QHN63861.1 is Class I, not Class IV — again matching the
config's prior `class_i=true` flag on that accession.

### Verdict on Task 1

**Confirmed, not contradicted.** Both flagged accessions (BAW35424.1, QHN63861.1) are genuinely
Class I in the primary literature, and Class I vs. Class IV is a real, independently-verifiable
domain-architecture distinction (intact CBM18 + loopful catalytic domain vs. truncated CBM18 +
loopless catalytic domain), not a labeling nitpick.

---

## Task 2: Mixing classes in one alignment is a real methodological problem

Standard phylogenetic/orthology practice: multiple sequence alignment and ancestral-state
reconstruction assume the aligned region is homologous across all sequences. When sequences differ
in domain architecture (here: presence/absence/truncation of an N-terminal CBM18 domain, plus
"loopful" vs. "loopless" catalytic-domain indels), the aligner is forced to either (a) introduce
large gap blocks that are not orthologous absences but domain-architecture differences, or (b)
misalign N-terminal catalytic-domain residues against CBM18 residues from the other class. Either
outcome corrupts column homology at exactly the positions ancestral reconstruction and convergence
counting depend on. This is exactly the scenario CLAUDE.md's own Phase 1 test gate is designed to
catch (tree topology / orthogroup composition checks), and Fukushima et al. 2017's convergence
method (ancestral reconstruction + counting derived-state matches on carnivorous-lineage branches)
requires a single consistent alignment column per site across the tree — a mixed Class I/IV
alignment would inject spurious "convergent" or "divergent" calls driven by domain architecture,
not by the actual catalytic-domain substitutions Fig. 3a reports. The GH19 Engineering Database
source above (PMC8547705) independently treats class as a first-order stratifying variable for
exactly this reason (sequence diversity and substrate scope are analyzed per class, not pooled).

---

## Gate decision: PASS

Both retrieved sources confirm the Class I/IV distinction as described in the task, and both
flagged accessions are genuinely Class I. Proceeded with the restructuring edit below.

---

## Edit applied to `config/enzyme_families.yaml`

### Before

Single tier1 entry `chitinases_gh19` contained Nepenthes/Cephalotus/Drosera_rotundifolia/
Dionaea/Sarracenia/outgroup Class IV-type accessions **and** two Class I accessions
(Drosera_adelae BAW35424.1, Drosera_binata QHN63861.1), plus a Drosera_capensis entry that used
a verified Class I accession (BDB33378.1) framed as a "proxy" for the missing Class IV gene.

### After

- **`chitinases_gh19_class_iv`** (tier1): everything previously in `chitinases_gh19` *except* the
  two Class I entries. `convergent_source` updated to explicitly state this is the Fukushima et
  al. 2017 Fig. 3a target class. `Drosera_capensis` now carries an honest `TODO` (see below)
  instead of the disguised Class I proxy. All other species blocks (Nepenthes_gracilis,
  Nepenthes_alata, Cephalotus_follicularis, Drosera_rotundifolia, Dionaea_muscipula,
  Arabidopsis_thaliana, Solanum_lycopersicum, Sarracenia_purpurea, Sarracenia_alata,
  Darlingtonia_californica, Heliamphora_ciliata, Utricularia_gibba, Pinguicula_moranensis,
  Byblis_filifolia, Drosera_spatulata, Nepenthes_mirabilis, Secale_cereale, and the Aldrovanda
  "REMOVED" note) preserved verbatim, unchanged.

- **`chitinases_gh19_class_i`** (tier1, new): Drosera_adelae (BAW35424.1) and Drosera_binata
  (QHN63861.1), each explicitly labeled "structural comparison only — not used for cross-lineage
  convergence detection." Drosera_capensis's BDB33378.1 (verified Class I) also moved here from
  the old proxy slot, since it is genuinely Class I and useful for structural comparison — it is
  no longer disguised as a Class IV substitute.

- **Drosera_capensis Class IV gap**: the old "use BDB33378.1 (class I) as proxy if class IV
  unavailable" framing (which also parenthetically mentioned an unused AZZ09188.1 alternate) is
  removed. Replaced with an honest `TODO` documenting that DCAP_0533 (the actual D. capensis
  Class IV gene identifier) exists only in an unsubmitted draft genome and has no public
  accession, following the same honesty convention as the existing Aldrovanda "REMOVED" notes.

No accession, TODO, or comment content was deleted — only reorganized/relabeled per the above.

### Snakefile / rule impact check

`Snakefile` lines 44–49 derive `FAMILIES = list(data.get("tier1", {}).keys())` generically at
parse time — confirmed by direct read. The two new keys (`chitinases_gh19_class_iv`,
`chitinases_gh19_class_i`) appear in `FAMILIES` automatically; no Snakefile change needed.
Verified this resolves correctly:

```
tier1 keys: ['chitinases_gh19_class_iv', 'chitinases_gh19_class_i', 'purple_acid_phosphatase',
             'rnase_t2', 'nepenthesins', 'aspartic_proteases_a1b_homology']
```

Grepped `workflow/rules/*.smk` and `workflow/scripts/*.py` for the literal string
`"chitinases_gh19"`:
- **`workflow/rules/*.smk`**: no matches.
- **`workflow/scripts/*.py`**: one match, `detect_convergence.py:410` — a `click.option` help
  string ("Enzyme family identifier (e.g. chitinases_gh19). Written into output TSV.") used only
  as an illustrative example in CLI help text, not a functional lookup. Does not break; harmless
  but now a slightly stale example (family names are passed as CLI args at runtime, not hardcoded).

**Additional dangling references found outside the requested grep scope** (config files, not
`.smk`/`.py`), which *would* have silently broken substrate lookups for the Class IV family after
the split, so fixed as part of this task:
- `config/config.yaml` — `fep.substrates_for_fep.chitinases_gh19: "GlcNAc4"` → renamed key to
  `chitinases_gh19_class_iv`.
- `config/substrates.yaml` — `GlcNAc4.for_family: chitinases_gh19` → updated to
  `chitinases_gh19_class_iv`.

**Found but left unfixed (out of scope, pre-existing, not part of the Snakemake pipeline)**:
`workflow/scripts/braker2/04_run_hmmer_scan.sh` is a standalone manual genome-mining helper
(uncommitted from a prior task) with a `FAMILY_HMMS`/`MIN_LEN` associative array keyed on
`"chitinases_gh19"` (and still `"nepenthesins"`, already stale from the prior A1B split in
audit/03). It does not read `enzyme_families.yaml` and is not invoked by `Snakefile` or any
`.smk` rule, so nothing breaks programmatically — but its labels/output filenames
(`*_chitinases_gh19.faa`) will need updating to `*_chitinases_gh19_class_iv.faa` (and a
`nepenthesins`/`aspartic_proteases_a1b_homology` split) by whoever next touches that script.
Flagged for the pipeline maintainer, matching the precedent in audit/03 for out-of-scope
dangling references.

---

## Verification

Independent re-check performed by a separate checking agent (did not make the original edit).
Did not trust the write-up above — re-fetched sources and re-read the current repo state fresh.

### 1. Independent literature re-retrieval

- **Yoneda et al. 2025** (`febs.onlinelibrary.wiley.com/doi/full/10.1002/2211-5463.70110`),
  fetched directly: confirms *D. adelae* chitinase **BAW35424.1** "belongs to class I of the
  GH19 family," defined by "the presence of a chitin-binding domain, also known as the hevein
  domain, at their N-terminal region." Crystal structure PDB **9JTR** (wild-type, 1.73 Å) plus
  9JTP (E167Q mutant, 1.57 Å) — matches the audit's citation exactly (only the catalytic domain
  was resolved experimentally; the hevein domain itself was AlphaFold2-modeled, a nuance not
  material to the class-I claim).
- **PMC8547705** ("The GH19 Engineering Database"), fetched directly (independent of the
  original WebSearch paraphrase): "Classes I and IV are linked to one accessory N-terminal
  carbohydrate binding module (CBM), whereas class II GH19s are characterized by the absence of
  a CBM"; "The sequences of class IV enzymes are shorter than classes I and II, resulting in a
  smaller number of subsites in the catalytic cleft"; loop architecture is split into "loopful"
  (class I/II) vs. "loopless" (class IV) categories. This is a *refinement*, not a contradiction,
  of the audit's "CBM18 with a deletion" framing: class IV retains a CBM (unlike class II, which
  has none) but a truncated/shortened one — exactly what "deletion" implies, not full absence.
  The audit's phrasing is accurate.
- A second independent WebSearch (not reusing the original agent's query) returned the same
  class I/II/IV definitional quote verbatim, corroborating both the domain-architecture
  distinction and the "loopful"/"loopless" terminology.
- Did not re-fetch the Planta 2025 (PMC11725546) *D. binata* source separately, per the task's
  "don't need to re-derive every detail" allowance, given the Yoneda + PMC8547705 results were
  unambiguous and consistent with the config's claims.

**Conclusion:** the Class I vs. Class IV domain-architecture distinction, and both accession
classifications (BAW35424.1, QHN63861.1 = Class I), are independently confirmed.

### 2. Fresh read of `config/enzyme_families.yaml`

Parsed the live file with PyYAML (not just eyeballed) after a full manual read:

```
tier1 keys: ['chitinases_gh19_class_iv', 'chitinases_gh19_class_i', 'purple_acid_phosphatase',
             'rnase_t2', 'nepenthesins', 'aspartic_proteases_a1b_homology']
class_iv species: ['Nepenthes_gracilis', 'Nepenthes_alata', 'Cephalotus_follicularis',
  'Drosera_rotundifolia', 'Drosera_capensis', 'Dionaea_muscipula', 'Arabidopsis_thaliana',
  'Solanum_lycopersicum', 'Sarracenia_purpurea', 'Sarracenia_alata', 'Darlingtonia_californica',
  'Heliamphora_ciliata', 'Utricularia_gibba', 'Pinguicula_moranensis', 'Byblis_filifolia',
  'Drosera_spatulata', 'Nepenthes_mirabilis', 'Secale_cereale']
class_i species: ['Drosera_capensis', 'Drosera_adelae', 'Drosera_binata']
class_i Drosera_capensis: ['BDB33378.1']
class_iv Drosera_capensis: ['TODO']
```

- Both `chitinases_gh19_class_iv` and `chitinases_gh19_class_i` exist as separate top-level
  `tier1` entries. File parses as valid YAML.
- Every species/accession from the described old combined `chitinases_gh19` entry is accounted
  for in exactly one of the two new entries: 18 species (incl. all TODO placeholders and the
  Aldrovanda "REMOVED" note, verbatim) under `_class_iv`, and the 3 Class I entries
  (Drosera_capensis/BDB33378.1, Drosera_adelae/BAW35424.1, Drosera_binata/QHN63861.1) under
  `_class_i`. No duplicate accession appears in both entries; no accession dropped.
- `Drosera_capensis` under `_class_iv` is a `TODO` with prose explicitly stating "BDB33378.1 was
  previously listed here as a 'proxy' for the missing class IV gene — that framing was WRONG,"
  and that BDB33378.1 "is now listed under chitinases_gh19_class_i instead ... structural
  comparison only." This is honest — it does not imply BDB33378.1 covers Class IV anywhere.

### 3. Fresh read of `config/config.yaml` and `config/substrates.yaml`

- `config/config.yaml` line 119: `fep.substrates_for_fep` key is `chitinases_gh19_class_iv:
  "GlcNAc4"` (confirmed via direct read and via PyYAML parse: `list(c['fep']['substrates_for_fep'].keys())`
  includes `chitinases_gh19_class_iv`, not the old bare `chitinases_gh19`).
- `config/substrates.yaml` line 27: `GlcNAc4.for_family: chitinases_gh19_class_iv` (confirmed via
  direct read and PyYAML: `s['substrates']['GlcNAc4']['for_family'] == 'chitinases_gh19_class_iv'`).
- Both fixes are present and correct.

### 4. Repo-wide grep for the literal string `chitinases_gh19`

Ran `Grep -n "chitinases_gh19"` (not restricted to `.smk`/`.py` as the original audit's grep
was) across the whole repo. Matches, beyond the two new `_class_iv`/`_class_i` keys themselves:

- `audit/CHANGELOG.md`, `audit/03_merops_restructure_and_neprosin_rescope.md`,
  `audit/04_gh19_class_split.md` (this file) — historical/documentation references describing
  the old key by name. Expected and harmless; these are audit records, not live config/code.
- `workflow/scripts/braker2/04_run_hmmer_scan.sh` (lines 30, 42, 128) — the flagged standalone
  manual helper script. Confirmed not read by `Snakefile` or any `.smk` rule (it's a
  hand-invoked genome-mining script). Correctly left unfixed and correctly flagged.
- `workflow/scripts/detect_convergence.py:410` — `--family` CLI option's `help=` string uses
  `chitinases_gh19` purely as an illustrative example in documentation text, not a functional
  lookup key. Already identified in the original audit; independently re-confirmed by reading
  the surrounding code (lines 400-419) — it's inert.

No other functional/live-code or config reference to the bare `chitinases_gh19` key was found.

### 5. `Snakefile` FAMILIES derivation

Read `Snakefile` lines 38-49 directly:

```python
def _load_tier1_families(path: str) -> list[str]:
    with open(path, encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    return list(data.get("tier1", {}).keys())

FAMILIES = _load_tier1_families("config/enzyme_families.yaml")
```

This is a generic `tier1` key enumeration with no hardcoded family names. Combined with the
PyYAML parse above showing both `chitinases_gh19_class_iv` and `chitinases_gh19_class_i` as
live `tier1` keys, `FAMILIES` will include both at Snakemake parse time exactly as claimed.

### Verdict: **CONFIRMED**

All five claims hold up under independent re-verification: the Class I/IV literature distinction
is real and correctly sourced (re-fetched both primary sources directly, one additional
independent WebSearch as cross-check); the enzyme_families.yaml split is complete, non-lossy,
non-duplicated, structurally valid, and the Drosera_capensis TODO is honest; both dangling
references in config.yaml and substrates.yaml are fixed correctly; no other live reference to
the bare `chitinases_gh19` key exists outside the one correctly-flagged braker2 helper script and
one harmless CLI help-text example; and the Snakefile's generic tier1-key derivation will pick up
both new family keys automatically. No discrepancies found.
