# CarnivorEnzyme

A reproducible Snakemake pipeline for structural and functional characterization of convergently evolved digestive enzymes across all major carnivorous plant lineages.

---

![Nepenthes rajah pitchers](https://upload.wikimedia.org/wikipedia/commons/thumb/c/c8/Nepenthes_rajah_Schwarz.jpg/1024px-Nepenthes_rajah_Schwarz.jpg)

*Nepenthes rajah, Mount Kinabalu, Borneo. Pitcher fluid pH ~2.5; secretome includes aspartic proteases, chitinases, phosphatases, and glutamic peptidases with no phylogenetic continuity to functionally analogous enzymes in distantly related carnivorous lineages.*

---

## Background

Carnivory has evolved independently at minimum six times within the angiosperms — in *Nepenthes*, *Drosera*, *Dionaea*, and *Aldrovanda* (Caryophyllales); *Cephalotus* (Oxalidales); *Sarracenia* (Ericales); and *Utricularia* and *Pinguicula* (Lamiales). Each lineage independently co-opted pre-existing hydrolase gene families — chitinases, acid phosphatases, RNases, and proteases — for extracellular digestive function. Fukushima et al. (2017, *Nat. Ecol. Evol.*) established that this convergence extends below the phenotypic level: the same amino acid substitutions recur at homologous alignment positions across phylogenetically independent lineages at a rate significantly exceeding neutral expectation. The structural and functional consequences of these substitutions remain uncharacterized. No cross-lineage three-dimensional comparison of any carnivorous plant digestive enzyme family has been published.

This project addresses that gap by integrating modern structure prediction, thermodynamic free-energy perturbation, and evolutionary constraint inference across five enzyme families and all major carnivorous lineages.

---

## Target Enzyme Families

All five families have published convergent substitutions (Fukushima et al. 2017) and are expressed in the digestive secretome of multiple independent lineages.

| Family | Class | Digestive role | Reference structure |
| ------ | ----- | -------------- | ------------------- |
| GH19 chitinases | Glycoside hydrolase 19 | Chitin hydrolysis from arthropod cuticle | 4J0L (*Secale cereale*, 1.75 Å) |
| Nepenthesins / droserasins | A1B aspartic protease | Bulk proteolysis; optimum pH 2–3 | No coordinates deposited |
| Neprosins | G3 glutamic peptidase | Prolyl endopeptidase; Pro-Xaa cleavage | 7ZVA/B/C (*N.* × *ventrata*, 1.80–2.35 Å) |
| Purple acid phosphatases | Fe-Zn/Fe-Mn binuclear metalloenzyme | Organic phosphate liberation | 1RDP (*Glycine max*, 2.65 Å) |
| RNase T2 (S-like) | Ribonuclease T2 | RNA hydrolysis; phosphate scavenging | 1IYC (*Pyrus communis*, 1.50 Å) |

Neprosins are the only glutamic peptidases identified in plants (Tiew & Goh 2022). Their crystal structures revealed a catalytic dyad mechanistically distinct from both aspartic and serine protease classes, and their strict prolyl endopeptidase activity has generated clinical interest in celiac disease therapy, where the immunodominant gliadin epitopes are proline-rich and resistant to human gastrointestinal proteases.

---

![Cephalotus follicularis](https://upload.wikimedia.org/wikipedia/commons/thumb/7/72/Cephalotus_follicularis2.jpg/800px-Cephalotus_follicularis2.jpg)

*Cephalotus follicularis, southwestern Australia — the sole carnivorous species in Oxalidales, whose nearest relatives are oxalis and starfruit. Its digestive enzyme sequences are independently derived from every other carnivorous lineage, yet share convergent substitutions at multiple alignment positions with Nepenthes, Drosera, and Dionaea. The Fukushima 2017 genome paper used Cephalotus as the primary case study.*

---

## Analytical Design

The pipeline executes ten phases on a Snakemake DAG targeting a SLURM HPC environment (AlphaFold3 on A100/H100 GPUs). The analysis produces not just hypotheses but thermodynamically grounded answers.

**Phylogenomics (Phases 1–3).** Protein sequences for 20 carnivorous and 6 outgroup species are retrieved from NCBI and UniProt, aligned with MAFFT L-INS-i, and trimmed with trimAl. IQ-TREE 3 infers per-family maximum-likelihood gene trees with ModelFinder Plus model selection and the `-ancestral` flag for marginal ancestral state reconstruction. Convergent positions are defined as alignment sites where the same derived residue was independently fixed on ≥2 carnivorous lineage branches, evaluated against a Poisson null parameterized by branch-specific substitution rates (posterior probability threshold ≥0.8; cf. Fukushima & Pollock 2023).

**ESM-2 fitness triage (Phase 3A).** Before committing to expensive downstream computation, each convergent substitution is scored with ESM-2 masked marginal log-likelihood (Lin et al. 2023, *Science* 379:1123). The log-likelihood ratio (LLR) quantifies whether the derived amino acid is sequence-contextually favored. Sites with LLR > 0 — consistent with positive selection — are prioritized for FoldX, CpHMD, and FEP.

**Structure prediction (Phase 4).** AlphaFold3 (5 seeds per Tier 1 target) predicts structures for all unique sequences; Chai-1 serves as a fallback. Validation requires TM-score >0.90 against the neprosin crystal structure 7ZVA (1.80 Å). Residues with pLDDT <50 are masked in all subsequent steps.

**Ancestral structure reconstruction (Phase 4A).** The IQ-TREE marginal ancestral state reconstruction (.state file) is parsed to recover the MAP sequence of the carnivore MRCA for each family. AlphaFold3 predicts the ancestral structure. TM-align comparisons between the ancestral and each modern carnivorous ortholog — restricted to convergent positions — test whether those sites are structurally displaced relative to background, which would implicate them in the evolutionary acquisition of digestive function.

**Convergent substitution classification (Phase 5).** Each convergent substitution is characterized by two orthogonal methods. FoldX 5.1 computes folding free energy change (ΔΔG) using the revised 2025 force field with explicit pH-dependent protonation, run at pH 2.5, 3.5, and 5.0. EVcouplings independently fits a Potts model of residue–residue coevolution and scores each substitution's evolutionary fitness consequence (ΔΔE). Comparing ΔΔG against ΔΔE classifies each convergent site mechanistically:

```text
                        EVmutation ΔΔE
                     favorable  |  unfavorable
                   _____________|_____________
 FoldX     neutral |  FUNCTIONAL |  NEUTRAL    |
 ΔΔG               |  GAIN       |  DRIFT      |
           ________|_____________|_____________|
           destab. |  STABILITY- |  DELETERIOUS|
                   |  FUNCTION   |  (false pos)|
                   |  TRADEOFF   |             |
                   |_____________|_____________|
```

No published study has applied a combined FoldX × EVmutation analysis to convergent substitutions in any organism. ESM-2 LLR provides a third independent axis: sites that are favorable on all three axes represent the strongest candidates for adaptive convergence.

**Alchemical free-energy perturbation (Phase 5D).** Sites classified as functional_gain or stability_function_tradeoff, and where ESM-2 LLR > 0, are advanced to GROMACS + pmx alchemical FEP. The thermodynamic cycle computes ΔΔG_bind — the change in substrate binding free energy attributable to the convergent substitution — with MBAR uncertainty estimation across 13 λ-windows × 5 ns × 3 replicates. A negative ΔΔG_bind is direct evidence that the convergent amino acid enhances substrate binding. Note: pmx does not support proline residue mutations; affected sites are excluded and logged.

**Constant-pH MD (Phase 5E).** GROMACS native λ-dynamics CpHMD (phbuilder setup; Gapsys et al. 2022, *JCTC*) determines whether convergent substitutions shift the protonation equilibrium of nearby titratable residues (Asp, Glu) at pitcher fluid pH values. This addresses the adaptation of active-site chemistry to the highly acidic secretion environment — a question uniquely relevant to carnivorous plant enzymes.

**Downstream analyses (Phases 6–8).** AutoDock Vina docks five substrates against all predicted structures. PDB2PQR + APBS computes electrostatic potential surfaces at pH 2.5, 3.5, and 5.0, enabling correlation of surface charge distribution with pitcher fluid pH across species. Salmon quantifies expression of target genes in the three lineages with available paired RNA-seq and genome assemblies.

All results are aggregated into a SQLite atlas and exposed through a Streamlit browser interface.

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

External dependencies (not conda-installable): FoldX 5+ (academic license, [foldxsuite.crg.eu](https://foldxsuite.crg.eu)), ADFR Suite 1.0 (Scripps Research), TM-align, ORCA 6 (QM/MM only; academic license free, [orcaforum.kofo.mpg.de](https://orcaforum.kofo.mpg.de)).

---

## References

- Fukushima K et al. (2017) Genome of the pitcher plant *Cephalotus* reveals genetic changes associated with carnivory. *Nat. Ecol. Evol.* 1:0059.
- Fukushima K, Pollock DD. (2023) Detecting macroevolutionary genotype–phenotype associations using error-corrected rates of protein convergence. *Nat. Ecol. Evol.*
- Del Amo-Maestro L et al. (2022) Molecular and in vivo studies of a glutamate-class prolyl-endopeptidase for coeliac disease therapy. *Nat. Commun.* 13:4446.
- Abramson J et al. (2024) Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 630:493–500.
- Lin Z et al. (2023) Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science* 379:1123–1130.
- Hopf TA et al. (2017) Mutation effects predicted from sequence co-variation. *Nat. Biotechnol.* 35:128–135.
- Botte M et al. (2025) FoldX force field revision with pH dependency and expanded training set. *Bioinformatics* 41:btaf064.
- Gapsys V et al. (2022) GROMACS native constant-pH MD (λ-dynamics). *JCTC* 18:6320.
- Gapsys V, Michielssens S, Seeliger D, de Groot BL. (2015) pmx: Automated protein structure and topology generation for alchemical perturbations. *JCTC* 11:4494.
- Shirts MR, Chodera JD. (2008) Statistically optimal analysis of samples from multiple equilibrium states. *JCP* 129:124105.

## License

MIT
