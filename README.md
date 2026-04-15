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

The pipeline executes eight sequential phases on a Snakemake DAG targeting a SLURM HPC environment (AlphaFold3 on A100/H100 GPUs).

**Phylogenomics (Phases 1–3).** Protein sequences for 20 carnivorous and 6 outgroup species are retrieved from NCBI and UniProt, aligned with MAFFT L-INS-i, and trimmed with trimAl. IQ-TREE 3 infers per-family maximum-likelihood gene trees with ModelFinder Plus model selection and the `-ancestral` flag for marginal ancestral state reconstruction. Convergent positions are defined as alignment sites where the same derived residue was independently fixed on ≥2 carnivorous lineage branches, evaluated against a Poisson null parameterized by branch-specific substitution rates (posterior probability threshold ≥0.8; cf. Fukushima & Pollock 2023).

**Structure prediction (Phase 4).** AlphaFold3 (5 seeds per Tier 1 target) predicts structures for all unique sequences; Chai-1 serves as a fallback. Validation requires TM-score >0.90 against the neprosin crystal structure 7ZVA (1.80 Å) before downstream analysis proceeds. Residues with pLDDT <50 are masked in all subsequent steps.

**Mutation effect analysis (Phase 5).** Each convergent substitution is characterized by two orthogonal methods. FoldX 5.1 computes folding free energy change (ΔΔG) using the revised 2025 force field, which adds explicit pH-dependent protonation of Asp, Glu, Lys, Arg, and Tyr — run at pH 2.5, 3.5, and 5.0 to span the measured range of pitcher fluid pH across species. EVcouplings independently fits a Potts model of direct residue–residue coevolution to a deep plant-wide MSA and scores each substitution's evolutionary fitness consequence (ΔΔE). Comparing ΔΔG against ΔΔE classifies each convergent site mechanistically:

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

To our knowledge, no published study has applied a combined FoldX × EVmutation analysis to convergent substitutions in any organism. The two methods have been benchmarked against each other (Livesey & Marsh 2023) but not deployed as orthogonal axes for mechanistic classification. Their complementarity — FoldX captures thermodynamic stability (r = 0.71, Delgado et al. 2025); EVmutation captures evolutionary constraint from coevolutionary couplings (median r = 0.56 across 34 DMS assays, Hopf et al. 2017) — is precisely what makes the quadrant framework informative: disagreement between the axes is the signal, not noise. A substitution that is thermodynamically neutral but evolutionarily favored implicates functional rather than stability-based selection; one that is destabilizing yet evolutionarily favored suggests a stability–function tradeoff under positive selection; one that is unfavorable on both axes is likely a false positive in the convergence detection step.

**Downstream analyses (Phases 6–8).** AutoDock Vina docks five substrates (GlcNAc₄, gliadin octapeptide PQPQLPYP, ApA dinucleotide, p-nitrophenyl phosphate, hemoglobin octapeptide VHLTPEEK) against all predicted structures. PDB2PQR + APBS computes electrostatic potential surfaces at pH 2.5, 3.5, and 5.0, enabling correlation of surface charge distribution with measured pitcher fluid pH across species. Salmon quantifies expression of target genes in the three lineages with available paired RNA-seq and genome assemblies (*N. gracilis* PRJDB8591, *Cephalotus* PRJDB4470, *Dionaea* PRJEB12493).

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

External dependencies (not conda-installable): FoldX 5+ (academic license, [foldxsuite.crg.eu](https://foldxsuite.crg.eu)), ADFR Suite 1.0 (Scripps Research), TM-align.

---

## References

- Fukushima K et al. (2017) Genome of the pitcher plant *Cephalotus* reveals genetic changes associated with carnivory. *Nat. Ecol. Evol.* 1:0059.
- Fukushima K, Pollock DD. (2023) Detecting macroevolutionary genotype–phenotype associations using error-corrected rates of protein convergence. *Nat. Ecol. Evol.*
- Del Amo-Maestro L et al. (2022) Molecular and in vivo studies of a glutamate-class prolyl-endopeptidase for coeliac disease therapy. *Nat. Commun.* 13:4446.
- Abramson J et al. (2024) Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 630:493–500.
- Hopf TA et al. (2017) Mutation effects predicted from sequence co-variation. *Nat. Biotechnol.* 35:128–135.
- Botte M et al. (2025) FoldX force field revision with pH dependency and expanded training set. *Bioinformatics* 41:btaf064.

## License

MIT
