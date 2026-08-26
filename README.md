# CarnivorEnzyme

A reproducible Snakemake pipeline testing whether **convergent amino acid substitutions in carnivorous plant digestive enzymes track adaptation to pitcher fluid pH**, and whether the evolutionary trajectory of those substitutions is dominated by **stability-function tradeoff, functional gain, or neutral drift**.

> Current status: restructured (May 2026). See [PROJECT_RESTRUCTURE_2026-05-12.md](PROJECT_RESTRUCTURE_2026-05-12.md) for the question-driven reorganization, [ENZYME_EVOLUTION_PAPER_STRUCTURE.md](ENZYME_EVOLUTION_PAPER_STRUCTURE.md) for the published-paper template this mirrors, [STATUS_2026-05-12.md](STATUS_2026-05-12.md) for project state, and [TASKS.md](TASKS.md) for the active work plan.

---

![Nepenthes rajah pitchers](https://upload.wikimedia.org/wikipedia/commons/thumb/c/c8/Nepenthes_rajah_Schwarz.jpg/1024px-Nepenthes_rajah_Schwarz.jpg)

*Nepenthes rajah, Mount Kinabalu, Borneo. Pitcher fluid pH ~2.5; secretome includes aspartic proteases, chitinases, phosphatases, and glutamic peptidases with no phylogenetic continuity to functionally analogous enzymes in distantly related carnivorous lineages.*

---

## The Question

> **Across five independent origins of carnivory, do convergent amino acid substitutions in digestive enzymes track adaptation to pitcher fluid pH (1.9–6.0), and which evolutionary trajectory — functional gain, stability-function tradeoff, or neutral drift — characterizes each convergent site?**
>
> *(Corrected 2026-08-22: carnivory within Caryophyllales — Nepenthes, Drosera, Dionaea, Aldrovanda — is one shared origin, not four separate ones; see `config/species.yaml` and `audit/02_carnivory_origin_reassignment.md`.)*

This is the open question explicitly flagged by Albert/Fukushima 2026 (*Trends in Genetics*) as the structural-functional follow-up to Fukushima et al. 2017's sequence-level convergence detection, framed against the quantitative pitcher-pH phenotype documented by Pavlovič 2025 (*New Phytologist*) and Freund 2022 (*Plant Physiology*).

### The five testable hypotheses

| H | Claim | Test |
| - | ----- | ---- |
| H1 | Convergent substitutions are biased toward acid-stabilizing positions | FoldX ΔΔG at pH 2.5 vs pH 7 across convergent vs random sites |
| H2 | Convergent substitution pattern correlates with species pitcher pH | PGLS regression of mean ΔΔG per species against pitcher pH |
| H3 | Convergent sites show positive selection on carnivorous lineages | PAML branch-site test; HyPhy aBSREL |
| H4 | Stability-function tradeoff is the dominant trajectory | FoldX × evolutionary-constraint quadrant distribution |
| H5 (optional, wet lab) | The pre-carnivory ancestor enzyme is inactive at pH 2.5 while modern is active | kcat/Km of resurrected ancestor vs modern at pH 2.5 / 5 / 7 |

---

## Background

Carnivory has evolved independently at least five times within the angiosperms — once in Caryophyllales (a single shared origin ancestral to *Nepenthes*, *Drosera*, *Dionaea*, and *Aldrovanda*; Fleck & Jobson 2023, *Plants* 12(19):3356; Fleischmann et al. 2018, in Ellison & Adamec, eds., *Carnivorous Plants: Physiology, Ecology, and Evolution*, OUP, ch. 3); *Cephalotus* (Oxalidales); *Sarracenia*, *Darlingtonia*, and *Heliamphora* (Ericales); *Utricularia* and *Pinguicula* (Lamiales/Lentibulariaceae); and *Byblis* (Lamiales/Byblidaceae). Each lineage independently co-opted pre-existing hydrolase gene families — chitinases, acid phosphatases, RNases, and proteases — for extracellular digestive function.

Fukushima et al. (2017, *Nat. Ecol. Evol.*) established that this convergence extends below the phenotypic level: the same amino acid substitutions recur at homologous alignment positions across phylogenetically independent lineages at a rate significantly exceeding neutral expectation. The structural and functional consequences of these substitutions remain uncharacterized — a question explicitly named as open by the Albert/Fukushima 2026 *Trends in Genetics* review.

Pitcher fluid pH varies across species from ~1.9 (*N. rafflesiana*) to ~6.0 (*N. ampullaria*), providing a **quantitative continuous trait** for phylogenetic comparative analysis. CarnivorEnzyme tests whether convergent substitutions track this phenotype with the rigor expected of recent enzyme evolution papers — Iñiguez et al. 2022 (*Sci. Adv.*) on Solanaceae Rubisco, the Rubisco convergence study (bioRxiv 2025), and Wang/Stiffler et al. 2025 (*PNAS*) on adaptive convergent evolution.

---

## Core Analytical Contribution

The central novel analysis is a **two-axis classification of convergent substitutions** built from the consensus of two stability methods (FoldX + one ML) and one evolutionary method (EVmutation or ProSST):

```text
                  Evolutionary constraint (EVmutation / ProSST LLR)
                       favorable  |  unfavorable
                     _____________|_____________
 Structural   neutral |  FUNCTIONAL |  NEUTRAL    |
 stability            |  GAIN       |  DRIFT      |
 (FoldX ΔΔG)  ________|_____________|_____________|
              destab. |  STABILITY- |  DELETERIOUS|
                      |  FUNCTION   |  (false pos)|
                      |  TRADEOFF   |             |
                      |_____________|_____________|
```

No published study has applied this classification to convergent substitutions in any organism. The closest precedents:

- **Gerasimavicius, Livesey & Marsh 2023** (*Protein Science* doi:10.1002/pro.4688) — established the "stable-but-inactive" mutation class via stability-vs-DMS analysis; the conceptual foundation of the quadrant.
- **DETANGO** (Ding et al., bioRxiv Feb 2026) decomposes stability/function for general variants but is not applied to convergent evolution.
- **Chugunov et al. 2025** (*RSC Chem. Biol.* 6:975) experimentally validated that FoldX and evolutionary scoring capture orthogonal properties in one enzyme, but did not implement the quadrant.
- The **Rubisco convergence study** (bioRxiv 2025.10.08.681247) used FoldX on convergent plant enzyme sites — but one family only, FoldX-only, no evolutionary axis.

CarnivorEnzyme is the first to apply the stability × evolution classification to **convergent substitutions across four cross-lineage enzyme families and five independent lineage origins** (plus two single-lineage *Nepenthes* families retained as structure-prediction methods benchmarks), with quantitative correlation against pitcher fluid pH.

---

## Target Enzyme Families (Tier 1)

> **Restructured 2026-08-22** (see `audit/03_merops_restructure_and_neprosin_rescope.md` and
> `audit/04_gh19_class_split.md`), with GH19 Class I moved out of tier1 on 2026-08-24 (see
> `audit/10_config_gaps_and_consistency_test.md`). Four tier1 families are valid cross-lineage
> convergence targets; GH19 Class I, nepenthesins, and neprosins all live under
> `methods_benchmark:`, not `tier1:`, since each is confined to a single carnivory origin
> (GH19 Class I's three *Drosera* species all share the one Caryophyllales origin; nepenthesins
> and neprosins are *Nepenthes*-only) and a single-origin family cannot produce a
> `min_lineages: 2` cross-lineage convergence result. GH19 Class I is additionally kept out of
> the Class IV alignment/tree because it has a different domain architecture (intact CBM18
> chitin-binding domain vs. Class IV's truncated/deleted CBM18).

| Family | Class | Digestive role | Reference structure | Convergence target? |
| ------ | ----- | -------------- | ------------------- | -------------------- |
| GH19 chitinases, Class IV (`chitinases_gh19_class_iv`) | Glycoside hydrolase 19 | Chitin hydrolysis from arthropod cuticle | 4J0L (*Secale cereale*, 1.75 Å) | Yes — primary Fukushima 2017 Fig. 3a target |
| GH19 chitinases, Class I (`methods_benchmark.chitinases_gh19_class_i`) | Glycoside hydrolase 19 | Chitin hydrolysis (intact CBM18 chitin-binding domain) | PDB 9JTR (*D. adelae*, 1.73 Å) | No — single-origin structural comparison only |
| Aspartic proteases, A1B-like by homology (`aspartic_proteases_a1b_homology`: droserasin, dionain-AP, *Cephalotus* AP, *Sarracenia* AP; no MEROPS holotype) | A1B aspartic protease | Bulk proteolysis; optimum pH 2–3 | 1B5F (*Cynara cardunculus*, 2.0 Å) | Yes — spans Drosera+Dionaea, Cephalotus, Sarracenia origins |
| Purple acid phosphatases (`purple_acid_phosphatase`) | Fe-Zn/Fe-Mn binuclear metalloenzyme | Organic phosphate liberation | 1RDP (*Glycine max*, 2.65 Å) | Yes |
| RNase T2, S-like (`rnase_t2`) | Ribonuclease T2 | RNA hydrolysis; phosphate scavenging | 1IYC (*Pyrus communis*, 1.50 Å) | Yes |
| Nepenthesins (`methods_benchmark.nepenthesins`; MEROPS A01.040, *Nepenthes*-only) | A1B aspartic protease | Bulk proteolysis; optimum pH 2–3 | 1B5F (*Cynara cardunculus*, 2.0 Å) | No — single-origin methods benchmark |
| Neprosins (`methods_benchmark.neprosins`; MEROPS U74.001, "unknown catalytic type" — structurally, not taxonomically, related to family G3 per Oda & Wlodawer 2023) | Glutamic peptidase | Prolyl endopeptidase; Pro-Xaa cleavage | 7ZVA/B/C (*N.* × *ventrata*, 1.80–2.35 Å) | No — single-origin methods benchmark |

Neprosins are the only glutamic peptidases identified in plants (Ting et al. 2022, *Plant Physiology and Biochemistry* 183:23–35; corrected 2026-08-22 from a prior "Tiew & Goh 2022" misattribution — see `audit/03_merops_restructure_and_neprosin_rescope.md`). Their crystal structures revealed a catalytic dyad mechanistically distinct from aspartic and serine proteases, and their strict prolyl endopeptidase activity has generated clinical interest in celiac disease therapy (gluten contains proline-rich epitopes resistant to human gastrointestinal proteases).

---

![Cephalotus follicularis](https://upload.wikimedia.org/wikipedia/commons/thumb/7/72/Cephalotus_follicularis2.jpg/800px-Cephalotus_follicularis2.jpg)

*Cephalotus follicularis, southwestern Australia — the sole carnivorous species in Oxalidales, whose nearest relatives are oxalis and starfruit. Its digestive enzyme sequences are independently derived from every other carnivorous lineage, yet share convergent substitutions at multiple alignment positions with Nepenthes, Drosera, and Dionaea.*

---

## Analytical Design

The pipeline executes nine phases on a Snakemake DAG targeting a SLURM HPC environment (AlphaFold 3 on A100/H100 GPUs at the University of Missouri Hellbender cluster).

### Phylogenomics (Phases 1–2A)

Protein sequences for 20+ carnivorous and 5 outgroup species are retrieved from NCBI and UniProt, aligned with MAFFT L-INS-i, and trimmed with trimAl. **IQ-TREE 3** infers per-family maximum-likelihood gene trees with ModelFinder Plus model selection and the `-ancestral` flag for marginal ancestral state reconstruction.

**Convergent positions** are defined as alignment sites where the same derived residue was independently fixed on ≥ 2 carnivorous lineage branches, evaluated against a Poisson null parameterized by branch-specific substitution rates (posterior probability threshold ≥ 0.8; cf. Fukushima & Pollock 2023).

### Phase 2B — ACEP embedding-distance convergence

In addition to site-level counting, ortholog pairwise distance in protein-language-model embedding space tests whether carnivorous lineages converge functionally without site-level overlap (Stiffler et al. 2025, *PNAS* doi:10.1073/pnas.2418254122).

### Phase 3A — Structure-aware PLM scoring

**ProSST** (Li et al., *NeurIPS* 2024) and **SaProt** (Su et al., *ICLR* 2024) provide the structure-aware evolutionary axis of the quadrant analysis. Both ingest the predicted AF3 / Chai-1 structure and a sequence and output a per-position log-likelihood ratio. ProSST is currently the highest-ranked open model on the ProteinGym benchmark (first to exceed Spearman 0.5 on the full substitution benchmark). For deep-MSA families (Neff ≥ 500), **EVE** (Frazer et al. 2021) provides a complementary MSA-based axis.

### Phase 4 — Structure prediction

**AlphaFold 3** (5 seeds per Tier 1 target) predicts structures for all unique sequences; **Chai-1** serves as a fallback. Validation requires TM-score > 0.90 against the neprosin crystal structure 7ZVA (1.80 Å). Residues with pLDDT < 50 are masked in all subsequent steps.

### Phase 4A — Ancestral structure reconstruction

The IQ-TREE marginal ancestral state reconstruction (.state file) is parsed to recover the MAP sequence of the carnivore MRCA for each family. AlphaFold 3 predicts the ancestral structure. TM-align comparisons between the ancestral and each modern carnivorous ortholog — restricted to convergent positions — test whether those sites are structurally displaced relative to background, which would implicate them in the evolutionary acquisition of digestive function.

### Phase 5 — Convergent substitution classification (THE CORE NOVEL OUTPUT)

Each convergent substitution is characterized by **two orthogonal axes**:

- **Structural stability (X-axis):** FoldX 5.1 (Botte et al. 2025) computes folding free energy change (ΔΔG) using the revised 2025 force field with explicit pH-dependent protonation, run at pH 2.5, 3.5, and 5.0 — matching the range of pitcher fluid acidity across species.
- **Evolutionary constraint (Y-axis):** ProSST and SaProt LLR for each convergent substitution. For deep-MSA families, EVE/EVcouplings ΔΔE provides an MSA-based check.

The **FoldX × ProSST scatter plot per family** (5 panels in Fig. 5 of the manuscript) is the primary novel analytical figure. Each convergent site is classified mechanistically (functional_gain, stability_function_tradeoff, neutral_drift, deleterious).

### Phase 5D — Alchemical FEP binding free energy

Sites in the `functional_gain` or `stability_function_tradeoff` quadrant with ProSST LLR > 0 are advanced to **GROMACS + pmx alchemical FEP** (capped at 20 sites total project-wide). The thermodynamic cycle computes ΔΔG_bind — the change in substrate binding free energy attributable to the convergent substitution — with MBAR uncertainty estimation across 13 λ-windows × 5 ns × 3 replicates. A negative ΔΔG_bind is direct evidence that the convergent amino acid enhances substrate binding.

> Note: pmx does not support proline residue mutations; affected sites are excluded and logged.

### Phase 5E — Constant-pH MD

**GROMACS native λ-dynamics CpHMD** (phbuilder setup; Gapsys et al. 2022, *JCTC*) determines whether convergent substitutions shift the protonation equilibrium of nearby titratable residues (Asp, Glu, His, Lys, Cys, Tyr) at pitcher fluid pH values. This addresses the adaptation of active-site chemistry to the highly acidic secretion environment — a question uniquely relevant to carnivorous plant enzymes.

### Downstream analyses (Phases 6–8)

AutoDock Vina docks five substrates against all predicted structures. PDB2PQR + APBS computes electrostatic potential surfaces at pH 2.5, 3.5, and 5.0, enabling correlation of surface charge distribution with pitcher fluid pH across species. Salmon quantifies expression of target genes in three lineages with available paired RNA-seq and genome assemblies.

### Phase 9 — Integration

All results are aggregated into a **SQLite atlas** with tables for enzymes, structures, convergence, foldx, prosst, evcouplings, docking, electrostatics, fep, cphmd, ancestral, and expression. Ten publication figures are generated by [generate_figures.py](workflow/scripts/generate_figures.py). A Streamlit browser interface exposes the atlas for community use.

---

## Lineage Coverage

> **Corrected 2026-08-22** (see `config/species.yaml` header and `audit/02_carnivory_origin_reassignment.md`):
> carnivory within Caryophyllales is **one shared origin**, not one per genus. Fleck & Jobson
> (2023, *Plants* 12(19):3356) and Fleischmann et al. (2018, in Ellison & Adamec, eds.,
> *Carnivorous Plants: Physiology, Ecology, and Evolution*, OUP, ch. 3) both describe *Nepenthes*,
> *Drosera*, *Dionaea*, and *Aldrovanda* as descending from a single common carnivorous ancestor.
> This table's origin numbering now matches `config/species.yaml`'s `carnivory_origin` field
> exactly (5 origins total, not 7).

| Origin | Order | Species | Genome status |
| ------ | ----- | ------- | ------------- |
| 1 | Caryophyllales | *Nepenthes gracilis*, *N. mirabilis*, *Drosera capensis* (12x), *D. spatulata*, *D. regia*, *Dionaea muscipula* (4x), *Aldrovanda vesiculosa* | GCA_030504385.1; PRJEB86749 (PacBio HiFi, 2025); draft / dodecaploid (*Drosera*); tetraploid / partial (*Dionaea*/*Aldrovanda*) |
| 2 | Oxalidales | *Cephalotus follicularis* | GCA_001941015.1 |
| 3 | Ericales | *Sarracenia purpurea* | bioRxiv 2025.12.26.696377 (chromosome-scale, Dec 2025) |
| 3 | Ericales | *Darlingtonia californica*, *Heliamphora ciliata* | Tarnita 2023 (raw assembly; BRAKER2 annotation in progress) |
| 4 | Lamiales (Lentibulariaceae) | *Utricularia gibba*, *Pinguicula gigantea*, *P. moranensis* | GCA_002189035.1; bioRxiv 2025.04.05.646448 |
| 5 | Lamiales (Byblidaceae) | *Byblis filifolia* | no published genome |

**Non-carnivorous outgroups:** *Arabidopsis thaliana*, *Vitis vinifera*, *Solanum lycopersicum*, *Hordeum vulgare*, *Secale cereale*.

**Polyploidy.** Three of the most important carnivorous genomes are polyploid (*N. gracilis* 10x, *D. muscipula* 4x, *D. capensis* 12x). The pipeline now annotates each accession with subgenome of origin where multiple homeologs exist, addressing the Albert/Fukushima 2026 *Trends in Genetics* open question about which subgenome copies retain convergent substitutions.

---

## Quick Start

```bash
git clone https://github.com/piercetaylor/CarnivorEnzyme.git
cd CarnivorEnzyme
conda env create -f environment.yml
conda activate carnivorenzyme

cp .env.example .env          # add NCBI_EMAIL and optionally NCBI_API_KEY
snakemake -n                  # dry-run to verify DAG
snakemake --use-conda phase1 --cores 4
snakemake --use-conda phase2 --executor slurm --profile config/slurm/
```

**External dependencies (not conda-installable):**

- FoldX 5.1+ — academic license at [foldxsuite.crg.eu](https://foldxsuite.crg.eu)
- ADFR Suite 1.0 — Scripps Research
- TM-align — [zhanggroup.org/TM-align](https://zhanggroup.org/TM-align)
- ORCA 6 (QM/MM only; disabled by default) — academic license free at [orcaforum.kofo.mpg.de](https://orcaforum.kofo.mpg.de)
- AlphaFold 3 — local Hellbender installation; weights under CC-BY-NC-SA 4.0

---

## Documentation Map

| File | Purpose |
| ---- | ------- |
| **README.md** (this file) | Public-facing overview and scientific framing |
| **[STATUS_2026-05-12.md](STATUS_2026-05-12.md)** | Current project state, implementation status, pivot direction |
| **[TASKS.md](TASKS.md)** | Active task list, ranked by tier and priority |
| **[LITERATURE_REVIEW_2026-05-12.md](LITERATURE_REVIEW_2026-05-12.md)** | Strategic methods-stack literature review (May 2026) |
| **[CLAUDE.md](CLAUDE.md)** | Method justifications, tool versions, build order — for AI coding agents |
| **[AGENTS.md](AGENTS.md)** | Operational execution guide per phase — for AI coding agents |
| **[PROGRESS.md](PROGRESS.md)** | Historical April 2026 progress (superseded by STATUS_2026-05-12.md) |
| **[temporary_TODO_42526.md](temporary_TODO_42526.md)** | Manual data retrieval tasks (BRAKER2, ENA, JGI, CDD) |
| **[literature_review/](literature_review/)** | Comprehensive April 2026 bibliography (10 topic areas) |

---

## Key References

### Foundational

- Fukushima K et al. (2017) Genome of the pitcher plant *Cephalotus* reveals genetic changes associated with carnivory. *Nat. Ecol. Evol.* 1:0059. doi:10.1038/s41559-016-0059
- Fukushima K, Pollock DD (2023) Detecting macroevolutionary genotype–phenotype associations using error-corrected rates of protein convergence. *Nat. Ecol. Evol.* 7:1156. doi:10.1038/s41559-023-02087-9
- Albert VA, Avila-Robledillo L, Fleck S, et al. (2026) Complexity and innovation in carnivorous plant genomes. *Trends in Genetics*. doi:S0168-9525(26)00034-X
- Del Amo-Maestro L et al. (2022) Molecular and in vivo studies of a glutamate-class prolyl-endopeptidase for coeliac disease therapy. *Nat. Commun.* 13:4446. doi:10.1038/s41467-022-32215-1 [PDB 7ZVA-C]

### Methods (2024–2026)

- Abramson J et al. (2024) Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 630:493
- Li M et al. (2024) ProSST: Protein sequence-structure transformer with disentangled attention. *NeurIPS 2024*
- Su J et al. (2024) SaProt: Protein language modeling with structure-aware vocabulary. *ICLR 2024*
- Frazer J et al. (2021) Disease variant prediction with deep generative models of evolutionary data (EVE). *Nature* 599:91
- Stiffler MA et al. (2025) Language models reveal a complex sequence basis for adaptive convergent evolution of protein functions. *PNAS* doi:10.1073/pnas.2418254122
- Botte M et al. (2025) FoldX force field revision with pH dependency. *Bioinformatics* 41:btaf064
- Gapsys V et al. (2022) GROMACS native constant-pH MD (λ-dynamics). *JCTC* 18:6320
- Gapsys V et al. (2015) pmx: alchemical perturbations. *JCTC* 11:4494
- Hopf TA et al. (2017) Mutation effects from sequence co-variation (EVmutation/EVcouplings). *Nat. Biotechnol.* 35:128

### Direct precedents (must cite, must distinguish)

- Ding et al. (2026) DETANGO: Deconvolving mutation effects on protein stability and function. *bioRxiv* 2026.02.03.703560
- Chugunov AO et al. (2025) FoldX vs EVmutation experimental orthogonality in endoglucanase II. *RSC Chem. Biol.* 6:975. doi:10.1039/D5CB00013K
- Rubisco convergence study (2025) FoldX on convergent plant enzyme sites. *bioRxiv* 2025.10.08.681247

### New genome resources (2025–2026)

- Saul F et al. (2023) Subgenome dominance in *Nepenthes gracilis*. *Nat. Plants* 9:1213
- Goh HH et al. (2025) *Nepenthes mirabilis* PacBio HiFi genome. *PLOS ONE* doi:10.1371/journal.pone.0322885
- Albert/Lindqvist (2025) *Sarracenia purpurea* chromosome-scale genome. *bioRxiv* 2025.12.26.696377

Full bibliography: [LITERATURE_REVIEW_2026-05-12.md](LITERATURE_REVIEW_2026-05-12.md) and [literature_review/MASTER_LITERATURE_REGISTRY.md](literature_review/MASTER_LITERATURE_REGISTRY.md).

---

## License

MIT
