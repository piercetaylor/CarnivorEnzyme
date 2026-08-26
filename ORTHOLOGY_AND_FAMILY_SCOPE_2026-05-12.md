# Orthology, Homology, and the Neprosin Exception

> Date: 2026-05-12
> Triggered by the question: "Are the enzymes used in carnivory across all 9 sets even the same? or similar evolutionarily?"
> This note documents a structural correction to the project framing. Companion to [ENZYME_EVOLUTION_PAPER_STRUCTURE.md](ENZYME_EVOLUTION_PAPER_STRUCTURE.md) and [PROJECT_RESTRUCTURE_2026-05-12.md](PROJECT_RESTRUCTURE_2026-05-12.md).

---

## 1. The Issue

The project as originally framed treats **five enzyme families** (GH19 chitinases, A1B aspartic proteases, neprosins, PAPs, RNase T2) as Tier 1 cross-lineage convergence targets. But the convergent-substitution question only works for enzymes that are **homologous across lineages** — that is, paralogs of a common ancestral gene that each carnivorous lineage independently recruited into a digestive role, then accumulated lineage-specific substitutions.

Of the five Tier 1 families, **four are homologous across lineages and one is not.**

---

## 2. Family-by-Family Verdict

| Family | Universal in angiosperms? | Cross-lineage orthology | Project role |
|--------|----------------------------|--------------------------|--------------|
| GH19 chitinases (class IV) | Yes (defense gene family) | Each carnivorous lineage co-opted an independent class IV ortholog | ✅ Cross-lineage convergence target |
| Purple acid phosphatases (PAP) | Yes (Metallophos superfamily) | All lineages have ancestral PAPs; acidic-optimum paralogs recruited | ✅ Cross-lineage convergence target |
| RNase T2 (S-like) | Yes (RNS1/2/3 in *Arabidopsis*) | All lineages have S-like RNase orthologs | ✅ Cross-lineage convergence target |
| A1B aspartic proteases (phytepsins) | Yes (universal in plants) | Most carnivorous lineages confirmed; *Utricularia* / *Pinguicula* / *Byblis* uncertain | ✅ Cross-lineage convergence target (where data exists) |
| **G3 glutamic peptidases (neprosins)** | **NO — known only from *Nepenthes*** | **No orthologs in non-*Nepenthes* carnivorous lineages** | ⚠️ **Not a cross-lineage convergence target** |

### Why neprosins are different

Ting et al. 2022 (*Plant Physiology and Biochemistry* 183:23-35, doi:10.1016/j.plaphy.2022.04.027; [CORRECTED 2026-08-24 — previously misattributed as "Tiew & Goh 2022," see `audit/03_merops_restructure_and_neprosin_rescope.md`]) characterized neprosins as founding a new plant glutamic peptidase family, found exclusively in *Nepenthes*. The crystal structures (7ZVA, 7ZVB, 7ZVC, 7ZU8; Del Amo-Maestro 2022, *Nat. Commun.* 13:4446) are from *Nepenthes* × *ventrata*.

The original Fukushima 2017 (*Nat. Ecol. Evol.* 1:0059) convergence paper **did not analyze neprosins** for cross-lineage convergence — precisely because no orthologs exist in *Drosera*, *Dionaea*, *Cephalotus*, *Sarracenia*, *Utricularia*, *Pinguicula*, or *Byblis*. The Fukushima analysis covered chitinases, aspartic proteases, RNase T2, PAPs, and a few other widely-conserved families.

---

## 3. The Orthology Requirement (Tighter Than Family-Level)

Even within the homologous families, the cross-lineage convergence analysis requires **paralog-level orthology**, not just family-level homology. Examples:

- *Arabidopsis* has 24+ GH19 chitinases across classes I, II, III, IV. Carnivorous plants recruited **class IV** copies (small, acidic-optimum, secreted) — not random GH19s.
- PAPs include cytoplasmic and secreted paralogs. Carnivorous plants recruited secreted, acidic-pH-optimum paralogs.
- A1 aspartic proteases include intracellular phytepsins (vacuolar) and secreted forms. Nepenthesins, droserasins, and dionains are the secreted A1B subfamily — distinguished by the **plant-specific insert (PSI)** domain.

Fukushima 2017's actual finding: **carnivorous plants independently recruit the same paralog class** (e.g., class IV chitinases) **and that class then accumulates the same convergent substitutions** in independent lineages.

This means the project's [detect_convergence.py](workflow/scripts/detect_convergence.py) analysis only works if:

1. The per-family alignment includes only orthologous paralogs (e.g., class IV GH19 only — not class I/II/III mixed in)
2. The gene tree topology shows independent recruitment per lineage (each carnivorous-lineage class-IV-paralog is sister to non-carnivorous class-IV-paralog outgroups, not to other carnivorous-lineage paralogs)
3. The ancestral state reconstruction at the deepest common ancestor of carnivorous + non-carnivorous class IV paralogs gives a meaningful "pre-carnivory" state

This is why the test gates in [AGENTS.md](AGENTS.md) Phase 2 require manual tree inspection. The orthology curation must be elevated to a hard prerequisite for each family.

---

## 4. Lineage × Family Coverage Reality

The realistic accession grid (per [config/enzyme_families.yaml](config/enzyme_families.yaml) as of April 2026, after the audit):

|                          | GH19 | PAP | RNase T2 | A1B AP | Neprosins |
|--------------------------|------|-----|----------|--------|-----------|
| 1. *Nepenthes* (Caryo)   | ✓    | ✓   | ✓        | ✓      | ✓ |
| 2. *Drosera* (Caryo)     | ✓    | ✓   | ✓        | ✓      | ✗ |
| 3. *Dionaea*/*Aldrovanda* (Caryo) | ✓ | partial | ✓ | ✓ | ✗ |
| 4. *Cephalotus* (Oxalidales) | ✓ | ✓   | ✓        | ✓ (fixed: GAV73214.1) | ✗ |
| 5. *Sarracenia* (Ericales)   | ✓ | ✓   | TODO     | ✓      | ✗ |
| 5b. *Darlingtonia*/*Heliamphora* | TODO (BRAKER2 pending) | TODO | TODO | TODO | ✗ |
| 6. *Utricularia* (Lamiales)  | TODO (JGI) | TODO (JGI) | TODO (JGI) | partial | ✗ |
| 6b. *Pinguicula* (Lamiales)  | TODO (BRAKER2) | TODO | TODO | TODO | ✗ |
| 7. *Byblis* (Lamiales/Byblidaceae) | ✗ no genome | ✗ no genome | ✗ no genome | ✗ no genome | ✗ |

**Concrete count once manual TODOs resolve:**

- GH19: 8 of 9 origins with data (only *Byblis* missing)
- PAP: 8 of 9 origins with data
- RNase T2: 8 of 9 origins with data
- A1B AP: 7 of 9 origins (uncertainty in *Utricularia* and *Pinguicula*)
- Neprosins: **1 of 9 origins** — *Nepenthes* only

**The cross-lineage convergence question is testable on 4 families × ~6–8 origins per family = 25–30 lineage × family analyses.** That is still statistically substantial.

---

## 5. Implications for the Manuscript

### What neprosins should be in v1

Not a cross-lineage convergence target. Their legitimate roles in this project:

1. **Methods benchmark.** 7ZVA at 1.80 Å is the only experimental structure for any carnivorous plant digestive enzyme. Used to validate AF3 (TM-score gate ≥ 0.90), validate FoldX (Pearson r > 0.65), and validate docking (PQPQLPYP into 7ZVC active site).
2. **Within-*Nepenthes* case study.** Multiple *Nepenthes* species, each with measured pitcher pH — can test whether neprosin substitutions correlate with within-genus pitcher pH variation. This is a within-clade comparative analysis, not a cross-origin convergence claim.
3. **Functional-convergence-by-different-fold story.** Neprosins (G3 glutamic) and aspartic proteases (A1B) both serve prolyl-endopeptidase / general proteolysis roles in their respective lineages. This is **fold-level convergence**, methodologically distinct from amino-acid-substitution convergence — needs Foldseek structural comparison, not detect_convergence.py.

### Two manuscript options

**Option 1 — Single paper, 4 homologous families × 9 lineages.**

- Title (working): *"Convergent substitutions in carnivorous plant digestive enzymes track adaptation to acidic secretion environments"*
- Neprosins in supplementary methods only (AF3 / FoldX validation benchmark)
- Story: H1–H5 hypotheses, 4 families, PGLS against digestive_secretion_pH
- Cleaner statistical story. Avoids the neprosin-orthology problem entirely.

**Option 2 — Companion papers.**

- Paper A: same as Option 1 (4 families cross-lineage)
- Paper B: within-*Nepenthes* neprosin pH-adaptation case study, using 7ZVA crystal benchmark
- Smaller, more focused. Could be a *Plant Cell* or *New Phytologist* companion piece.

---

## 6. Documents to Update

| Document | Required change |
|----------|------------------|
| [README.md](README.md) | Move neprosins out of "Tier 1 cross-lineage convergence" into "methods benchmark / within-*Nepenthes* case" section |
| [CLAUDE.md](CLAUDE.md) Section 3 (Tier 1 families) | Re-tier neprosins as "methods benchmark" not cross-lineage convergence target |
| [config/enzyme_families.yaml](config/enzyme_families.yaml) | Add explicit comment on neprosin section: "Within-*Nepenthes* case + AF3/FoldX benchmark — NOT cross-lineage convergence" |
| [ENZYME_EVOLUTION_PAPER_STRUCTURE.md](ENZYME_EVOLUTION_PAPER_STRUCTURE.md) | Section 1 Q1: clarify that Fukushima 2017's actual analysis excluded neprosins because they are *Nepenthes*-specific |
| [PROJECT_RESTRUCTURE_2026-05-12.md](PROJECT_RESTRUCTURE_2026-05-12.md) Section 1 | Update H1–H5 to specify "across 4 families" not "across 5 families" |
| [STATUS_2026-05-12.md](STATUS_2026-05-12.md) Section 3 (Target Enzyme Families) | Add a "methods benchmark" tier with neprosins |
| [TASKS.md](TASKS.md) | T1.1 (Phase 2 execution) — add explicit orthology-curation step per family before convergence detection |

---

## 7. The Orthology Curation Step (Missing from Pipeline)

Currently the pipeline assumes that whatever sequences are in `results/sequences/{family}/combined.fa` are all orthologous within the family. **This is not enforced.** A class I GH19 from *Arabidopsis* and a class IV GH19 from *Nepenthes* could end up in the same alignment, and the convergence analysis would be polluted.

Required addition: **Section B0 — Orthology curation.**

Before alignment, for each family, the pipeline should:

1. Identify the relevant paralog class (e.g., class IV GH19; secreted A1B AP with PSI domain; cytoplasmic RNase T2 of the S-like clade; secreted-acidic PAP).
2. Either:
   - Filter combined.fa by HMM scan against the paralog-defining domain (e.g., PSI domain for A1B AP, chitin-binding signature for class IV GH19), OR
   - Run a preliminary tree on the broader family, identify the target clade, and prune to that clade
3. Verify clade composition before alignment (manual or automated test gate)

This step exists implicitly as a manual test gate in Phase 2 (see [PROGRESS.md](PROGRESS.md) test gates 2a–2c) but should be promoted to an automated script — `workflow/scripts/curate_orthologs.py` — with a clear input/output contract.

---

## 8. Sources

- Fukushima K et al. (2017) Genome of the pitcher plant *Cephalotus* reveals genetic changes associated with carnivory. *Nat. Ecol. Evol.* 1:0059
- Ting TY, Baharin A, Ramzi AB, Ng CL, Goh HH (2022) Neprosin belongs to a new family of glutamic peptidase based on in silico evidence. *Plant Physiol. Biochem.* 183:23-35. doi:10.1016/j.plaphy.2022.04.027 [CORRECTED 2026-08-24 — previously misattributed as "Tiew TY, Goh HH... 183:23"; see `audit/03_merops_restructure_and_neprosin_rescope.md`]
- Del Amo-Maestro L et al. (2022) Molecular and in vivo studies of a glutamate-class prolyl-endopeptidase for coeliac disease therapy. *Nat. Commun.* 13:4446 [PDB 7ZVA-C]
- Athauda SBP et al. (2004) Enzymic and structural characterization of nepenthesin. *Biochem. J.* 381:295 [A1B AP characterization]
- Butts CT, Bierma JC, Martin RW (2016) Novel proteases from the genome of *Drosera capensis*. *Proteins* 84:1517 [A1B in *Drosera*]
- Ting et al. 2022 — explicitly states neprosin family is unique to *Nepenthes*
