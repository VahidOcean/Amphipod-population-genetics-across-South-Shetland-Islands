# Population genetics of *Cheirimedon femoratus* across the South Shetland Islands

## Effective population size (Ne) estimation pipeline

This repository contains scripts for estimating effective population size (Ne) using the linkage disequilibrium (LD) method from GBS/RADseq SNP data.

**Species:** *Cheirimedon femoratus* (amphipod)  
**Sampling locations:** Four populations across the South Shetland Islands, Antarctica  
&emsp;— DICF (Deception Island), KGCF (King George Island), LICF (Livingston Island), SICF (Snow Island)

---

## Repository structure

```
.
├── 01_prepare_data.sh          # VCF → per-population PLINK PED/MAP files
├── Ne_estimation_LDmethod.R    # LD-based Ne estimation (Waples 2006)
└── README.md
```

---

## Requirements

| Software | Version | Link |
|---|---|---|
| PLINK | 1.9 | https://www.cog-genomics.org/plink/ |
| R | ≥ 4.0 | https://www.r-project.org/ |
| R: adegenet | ≥ 2.1 | `install.packages("adegenet")` |
| R: ggplot2 | ≥ 3.4 | `install.packages("ggplot2")` |
| R: dplyr | ≥ 1.1 | `install.packages("dplyr")` |

Install R packages:
```r
install.packages(c("adegenet", "ggplot2", "dplyr"))
```

---

## Input files

### VCF file
Standard VCF format. Produced by Stacks `populations` module or equivalent.  
Sample IDs must match those in the population map.

### Population map (`Amphi.txt`)
Tab or space-separated, no header:
```
DICF_001    DICF
DICF_002    DICF
KGCF_001    KGCF
LICF_001    LICF
SICF_001    SICF
```

---

## Usage

### Step 1 — Prepare per-population PED/MAP files

```bash
# On NeSI HPC
module load PLINK/1.09b6.16-GCC-9.2.0

bash 01_prepare_data.sh \
  --vcf    final.recode.vcf \
  --popmap Amphi.txt \
  --outdir ne_data \
  --maf    0.05 \
  --geno   0.10 \
  --mind   0.20
```

**Output structure:**
```
ne_data/
  plink_all/
    all_samples.bed/.bim/.fam    (all individuals, unfiltered)
    all_filtered.bed/.bim/.fam   (quality filtered)
  populations/
    DICF/  DICF.ped  DICF.map
    KGCF/  KGCF.ped  KGCF.map
    LICF/  LICF.ped  LICF.map
    SICF/  SICF.ped  SICF.map
  logs/
```

### Step 2 — Run Ne estimation

```bash
# On NeSI HPC
module load R/4.3.1-gimkl-2022a

Rscript Ne_estimation_LDmethod.R \
  --datadir ne_data/populations \
  --outdir  ne_results
```

**Output files:**
```
ne_results/
  Ne_LDmethod_summary.csv       # Ne table: all populations × all Pcrit values
  Ne_LDmethod_plot.pdf          # Comparison plot (Pcrit = 0.05)
  Ne_LDmethod_sensitivity.pdf   # Sensitivity across Pcrit thresholds
```

---

## Methods summary

**Data preparation (`01_prepare_data.sh`)**  
VCF files were filtered using PLINK v1.9 with a minimum minor allele frequency of 0.05, a maximum of 10% missing genotypes per SNP, and a maximum of 20% missing data per individual. Per-population PED/MAP files were generated for Ne estimation.

**Ne estimation (`Ne_estimation_LDmethod.R`)**  
Contemporary Ne was estimated using the linkage disequilibrium method (Hill 1981; Waples 2006). Pairwise r² values were calculated across a random subsample of 2,000 SNPs per population (seed = 42). The Waples (2006) bias correction was applied:

&emsp;r²_drift = mean(r²) − 1/n − 3.19/n²

Ne was estimated as:

&emsp;Ne = 1 / (3 × r²_drift)

with 95% confidence intervals derived by block jackknife resampling over 200 SNP-pair blocks. Analyses were performed at three minimum allele frequency thresholds (Pcrit = 0.05, 0.02, 0.01) to assess sensitivity to rare alleles.

---

## Results

| Population | N | Ne (Pcrit 0.05) | 95% CI | Note |
|---|---|---|---|---|
| DICF | 22 | 169.8 | 151.9–187.8 | Reliable |
| KGCF | 19 | 798.7 | 340.1–1257.4 | Wide CI — small sample |
| LICF | 16 | 86.0 | 78.1–93.9 | Smallest Ne |
| SICF | 13 | 103.8 | 56.0–151.6 | Consistent across Pcrit |

---

## Notes on GONE (not used)

GONE (Santiago et al. 2020) was initially attempted for temporal Ne reconstruction but was found to be unsuitable for this dataset. GBS data typically produces 1–3 SNPs per RAD locus/scaffold (mean: 2.15 SNPs/scaffold across 2,968 scaffolds), and GONE requires multiple SNPs per chromosome to model LD decay across physical distance. LD-based Ne estimation (NeEstimator approach) is the appropriate method for fragmented GBS assemblies.

---

## Citation

If you use these scripts, please cite:

**Method:**  
Waples RS (2006) A bias correction for estimates of effective population size based on linkage disequilibrium at unlinked gene loci. *Conservation Genetics* 7:167–184.  
Hill WG (1981) Estimation of effective population size from data on linkage disequilibrium. *Genetical Research* 38:209–216.

**Software:**  
Chang CC et al. (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience* 4:7.  
Jombart T & Ahmed I (2011) adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. *Bioinformatics* 27:3070–3071.

---

## Contact

For questions about this pipeline, please open a GitHub issue.
