# Manuscript Framing Review — Has This Question Been Posed?

> Date: 2026-05-12
> Triggered by the question: "Has the proposed title/question even been posed in the literature, or is it relevant?"
> Companion to [ORTHOLOGY_AND_FAMILY_SCOPE_2026-05-12.md](ORTHOLOGY_AND_FAMILY_SCOPE_2026-05-12.md) and [PROJECT_RESTRUCTURE_2026-05-12.md](PROJECT_RESTRUCTURE_2026-05-12.md)

---

## The Proposed Question and Title

**Working title (May 2026 draft):**

> *"Convergent substitutions in carnivorous plant digestive enzymes track adaptation to acidic secretion environments"*

**Question:** Across 9 independent origins of carnivory, do convergent amino acid substitutions in 4 homologous digestive enzyme families (GH19 chitinases, A1B aspartic proteases, RNase T2, PAPs) track adaptation to digestive-secretion pH (1.9–6.0), and which evolutionary trajectory characterizes each substitution?

---

## 1. Has This Question Been Posed?

**Stated implicitly across multiple papers but never tested with phylogenetic comparative methods at multi-family scale.** Three layers of prior art exist.

### Layer A — "Carnivorous enzymes are acid-active" (settled, single-enzyme level)

| Paper | Finding |
|-------|---------|
| Athauda et al. 2004 (*Biochem. J.* 381:295, doi:10.1042/BJ20031575) | Nepenthesin I/II from *N. distillatoria* optimally active pH 2–3; disulfide-rich stabilization |
| Buch et al. 2013/2014 (*J. Biol. Chem.*) | Dionain (cysteine protease) optimally active pH 3.6 |
| Schräder et al. 2017 (*MCP* 16:1162) | Neprosin glutamate-class endopeptidase active pH 2–4 |
| Hatano & Hamada 2012 (*J. Proteome Res.*) | *N. alata* fluid proteome with low-pH kinetics |
| Ravee, Salleh & Goh 2018 (*PeerJ* 6:e4914, doi:10.7717/peerj.4914) | Review of digestive enzymes by family, no phylogenetic correlation |

**At the single-enzyme, single-lineage level, acid activity is published fact.** This is what a reviewer means when they say "we already knew that."

### Layer B — "Convergent substitutions exist" (established, Fukushima 2017)

| Paper | Finding |
|-------|---------|
| Fukushima et al. 2017 (*Nat. Ecol. Evol.* 1:0059, doi:10.1038/s41559-016-0059) | Identified convergent AA substitutions in chitinases, PAPs, RNase T2 across *Cephalotus*, *Nepenthes*, *Drosera*, *Dionaea*. Mapped on homology models. **Did not test pH correlation.** **Did not apply thermodynamic or evolutionary-constraint quantification.** |

### Layer C — "These two facts should be linked" (proposed as open question, never formally answered)

| Paper | Statement |
|-------|-----------|
| Pavlovič 2025 (*New Phytol.* 247:2581, doi:10.1111/nph.70229) | Explicitly names structural/functional consequences of convergent substitutions as unresolved |
| Albert/Fukushima 2025 NSF review (PAR 10636585) | Same |
| Adlassnig, Peroutka, Lendl 2011 (*Ann. Bot.* 107:181, doi:10.1093/aob/mcq238) | Cross-lineage pitcher pH compilation (1.9–6.0). Provides the trait values needed for a PGLS test, but **does not perform the test.** |

**Verdict on Layer C:** The pH-tracking hypothesis is stated implicitly but has **never been formally tested with phylogenetic comparative methods across multiple families.** Pavlovič 2025 is the closest to a published hypothesis statement; it does not propose a methodological route to test it.

---

## 2. Is the Question Still Interesting?

**Partially. There is a genuine, defensible scientific opening — but the proposed title risks framing it as a derivative claim.**

### The reviewer concern that will appear

> "*Athauda 2004 already showed nepenthesins are acid-active. Schräder 2017 showed neprosins are acid-active. Buch 2014 showed dionain is acid-active. We have known for two decades that carnivorous plant enzymes operate at low pH. What new insight does the convergence-vs-pH framing provide?*"

This is the single strongest objection a reviewer at *Nature Ecology & Evolution* or *eLife* will raise.

### The legitimate counter-argument

The proposed analysis answers a question the prior literature did **not** ask:

> "**Which specific substitutions** implement that acid adaptation, **is the same mechanism reused across 9 independent origins** of carnivory, and **what fraction** of convergent substitutions in carnivorous plant enzymes are explained by pH adaptation vs by substrate-recognition shifts vs by stability-function tradeoffs vs by neutral drift?"

This is a **mechanism-of-convergence question**, not an **existence-of-acid-activity question**. The distinction matters and must be foregrounded in the abstract.

### What's actually the hot question in 2026

The Pavlovič 2025 / Albert-Fukushima 2025 / Saul 2023 *Nat. Plants* triad has set the field's attention on:

| Hot topic (2025–2026) | Citation anchor |
|-----------------------|------------------|
| Jasmonate signaling co-option | Pavlovič 2025 (*New Phytol.* 247:2581) |
| Host–microbe interactions in pitcher fluid | Hayashi et al. 2023 (*Appl. Environ. Microbiol.* doi:10.1128/aem.00812-23, N-fixing bacteria in less-acidic pitchers) |
| Gene duplication, polyploidy, subgenome dominance | Saul et al. 2023 (*Nat. Plants* 9:1213) on *N. gracilis* 10x |
| Substrate-recognition (neprosin gluten cleavage applications) | Ting et al. 2026 (*J. Future Foods*); Wall et al. 2026 (*ChemBioChem*) |

**Pure pH adaptation is NOT the hot 2026 question.** The proposed title positions the work against a question the field is not currently asking about.

---

## 3. The Closest Comparable Published Study

**Rubisco convergence/thermal adaptation** (bioRxiv 2025.10.08.681247) is the **methodological template**:

- Convergent substitutions across 4 plant clades
- FoldX RepairPDB + BuildModel
- PGLS correlation with growing-season temperature
- Single family (RbcL)
- FoldX-only stability axis

**No equivalent exists for digestive enzymes vs pH.** CarnivorEnzyme as currently designed would be the digestive-enzyme analogue, with three meaningful methodological advances:

1. **Multi-family scope** (4 enzyme families vs 1)
2. **Dual-axis mutation effect classification** (FoldX × evolutionary constraint — neither the Rubisco paper nor any carnivore-enzyme paper has done this)
3. **Modern structure prediction** (AF3/Chai-1 vs homology modeling)

---

## 4. Reviewer Predictions

### MBE / GBE reviewer (sympathetic)

> "The PGLS + convergence framework on a 9-lineage system is methodologically sound. The FoldX × EVmutation quadrant approach is genuinely novel. Concern: with only ~9 species and large within-species variance in pitcher pH (Adlassnig 2011), statistical power to detect pH-trait correlation is weak. Authors should report power analysis and treat the pH-PGLS as exploratory."

**Likely outcome:** accept with revisions. Probability ~60%.

### NEE reviewer (skeptical)

> "The acid-active phenotype of carnivorous plant enzymes is established (Athauda 2004; Schräder 2017). The paper needs to articulate clearly what is learned **beyond** confirming this at scale. The convergence-vs-pH framing risks being circular — these enzymes operate in acidic fluid, so finding 'convergent substitutions track pH' is unsurprising. The interesting result would be the **opposite**: which convergent sites are **not** explained by pH (i.e., substrate or stability)?"

**Likely outcome:** revise-and-resubmit, possible reject. Probability of acceptance without reframing: ~20–30%.

### eLife reviewer

Same NEE-style critique. Plus an additional concern: "Where is the wet-lab validation? Modern enzyme evolution at eLife typically includes at least one resurrected ancestor with kinetic characterization."

**Likely outcome:** desk reject probable without wet lab.

---

## 5. The Selling Angle That Works in 2026

**Reframe the question from "track pH adaptation" to "trajectory atlas with pH as one diagnostic axis."**

### Proposed revised title

> *"Stability, evolution, and pH adaptation: a structural atlas of convergent substitutions in carnivorous plant digestive enzymes"*

Or, more provocatively:

> *"How much of convergent enzyme evolution in carnivorous plants is explained by pH adaptation?"*

### Why this framing wins

- **Centers the quadrant analysis** (the genuinely novel methodological contribution) instead of pH (the already-known phenotype).
- **Frames pH as a hypothesis to test**, not a foregone conclusion. The "*how much*" framing accommodates the possibility that pH explains only 30–40% of the convergent signal, with the rest going to substrate specificity, stability, or drift — which is the **more interesting** scientific story.
- **Acknowledges Athauda 2004 / Schräder 2017** as the prior art (Methods + Discussion) instead of competing with them.
- **Pre-empts the "we already knew that" objection** by making the *fraction* of pH-explained convergence the headline finding, whether high or low.

### Restructured H1–H5

The hypotheses remain the same, but their narrative role changes:

| H | Claim | Role in revised framing |
|---|-------|--------------------------|
| H1 | Convergent substitutions are biased toward acid-stabilizing positions at pH 2.5 | Primary test, **fraction quantified** |
| H2 | Convergent substitution pattern correlates with species digestive-secretion pH (PGLS) | Diagnostic for pH-driven adaptation |
| H3 | Convergent sites show positive selection on carnivorous lineages | Tests whether convergence is adaptive vs neutral |
| H4 | Stability-function tradeoff is the dominant trajectory | **Core trajectory atlas — primary result** |
| H5 (optional) | Ancestor enzyme inactive at pH 2.5; modern active | Mechanism validation |

The change is in emphasis, not content. **H4 becomes the headline; H1+H2 become the diagnostic test of one trajectory within H4.**

---

## 6. The Unique Gap That CarnivorEnzyme Actually Fills

A truthful gap statement that survives reviewer scrutiny:

> "No prior work combines (a) cross-lineage scope across ≥4 independent origins of carnivory, (b) multi-family scope across ≥4 homologous digestive enzyme families, (c) modern *ab initio* structure prediction (AF3) replacing homology models, (d) dual-axis mutation-effect classification combining structural stability (FoldX) and evolutionary constraint (EVmutation / ProSST), and (e) phylogenetic comparative testing against environmental traits including pH. The Rubisco temperature-adaptation study has (a)+(c)+(e) for one enzyme; Fukushima 2017 has (a)+(b) with homology models only; Butts et al. 2016 has (c) for one family in one lineage. **The methodological combination is novel; the individual components are not.**"

This statement is defensible. The pure pH-adaptation framing is not.

---

## 7. Risk Acknowledgments (To Discuss in Manuscript)

These caveats should be explicit in Discussion to pre-empt review:

| Caveat | Mitigation |
|--------|-----------|
| Only 9 lineages — weak PGLS statistical power | Report 95% CI; treat pH-PGLS as exploratory; use bootstrap |
| Within-species pH variance (Adlassnig 2011) | Use species median pH; report sensitivity to pH choice |
| Acid activity already known for individual enzymes | Cite Athauda 2004 / Schräder 2017 explicitly in Intro; position contribution as "quantitative trajectory atlas" not "acid activity discovery" |
| Some lineages have microbial-dominated digestion (*Sarracenia*, *N. ampullaria*) | Add `digestion_mode` covariate to PGLS; or restrict to plant-secreted enzyme producers |
| No wet-lab validation in v1 | Target plant-friendly journals (Plant Cell, New Phytologist) or commit to Outcome B with one wet-lab experiment |
| Computational ΔΔG carries 1–2 kcal/mol error | Report distributions and rankings, not absolute energies |
| Predicted structures (AF3) not crystals | Cite Akdel 2022 NSMB; enforce pLDDT ≥ 70 |

---

## 8. Concrete Recommendations

1. **Reframe the title.** Move the headline from "track pH adaptation" to "stability-function trajectory atlas with pH as one diagnostic." Two candidate titles:
   - *"Stability, evolution, and pH adaptation: a structural atlas of convergent substitutions in carnivorous plant digestive enzymes"*
   - *"How much of convergent enzyme evolution in carnivorous plants is explained by pH adaptation? A four-family, nine-origin trajectory atlas."*

2. **Update [README.md](README.md), [STATUS_2026-05-12.md](STATUS_2026-05-12.md), [PROJECT_RESTRUCTURE_2026-05-12.md](PROJECT_RESTRUCTURE_2026-05-12.md) §1** to reflect the revised framing. Keep H1–H5 unchanged; change the narrative role of each.

3. **Cite Athauda 2004, Schräder 2017, Buch 2014 explicitly in the Introduction** as prior art that established acid activity. State that CarnivorEnzyme builds on this by asking the *fraction-of-convergence-explained-by-pH* question, not the *do-carnivorous-enzymes-work-at-low-pH* question.

4. **Position against the Rubisco convergence paper** (bioRxiv 2025.10.08.681247) as the methodological template, while highlighting the four-family and dual-axis extensions.

5. **Acknowledge in Discussion** that the field's hot questions in 2026 are jasmonate signaling, microbial mutualism, and polyploidy/subgenome dominance — and that CarnivorEnzyme provides the structural-thermodynamic foundation that those questions will eventually need.

6. **Realistic journal target:** **MBE** or **GBE** as the default. *Nature Ecology & Evolution* only with one wet-lab figure (Outcome B). *eLife* only with wet lab AND substantial reframing.

---

## 9. Sources

- Athauda SBP et al. (2004) Enzymic and structural characterization of nepenthesin. *Biochem. J.* 381:295. doi:10.1042/BJ20031575
- Buch F et al. (2013/2014) Dionain characterization. *J. Biol. Chem.* (multiple)
- Schräder CU et al. (2017) Neprosin selective prolyl endoprotease. *Mol. Cell. Proteomics* 16:1162
- Hatano N, Hamada T (2012) Proteome analysis of pitcher fluid of *Nepenthes*. *J. Proteome Res.*
- Ravee R, Salleh F, Goh HH (2018) Discovery of digestive enzymes in carnivorous plants. *PeerJ* 6:e4914
- Fukushima K et al. (2017) *Cephalotus* genome. *Nat. Ecol. Evol.* 1:0059
- Pavlovič A (2025) Diversity in digestion in carnivorous plants. *New Phytol.* 247:2581. doi:10.1111/nph.70229
- Albert/Fukushima 2025 NSF review (PAR 10636585) — "Complexity and Innovation in Carnivorous Plant Genomes"
- Adlassnig W, Peroutka M, Lendl T (2011) Traps of carnivorous pitcher plants as a habitat. *Ann. Bot.* 107:181. doi:10.1093/aob/mcq238
- Saul F et al. (2023) Subgenome dominance in *Nepenthes gracilis*. *Nat. Plants* 9:1213
- Hayashi N et al. (2023) N-fixing bacteria in less-acidic Nepenthes. *Appl. Environ. Microbiol.* doi:10.1128/aem.00812-23
- Ting TY et al. (2026) Bioengineered yeast with neprosin for gluten detoxification. *J. Future Foods* 6:1227
- Rubisco convergence/thermal adaptation (2025) bioRxiv 2025.10.08.681247
- Butts CT, Bierma JC, Martin RW (2016) Novel proteases from *Drosera capensis*. *Proteins* 84:1517. doi:10.1002/prot.25095
- Akdel M et al. (2022) AlphaFold2 community assessment. *Nat. Struct. Mol. Biol.* 29:1056. doi:10.1038/s41594-022-00849-w
