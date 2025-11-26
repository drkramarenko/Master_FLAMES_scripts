# Leveraging shared and opposing genetic mechanisms in heritable cardiomyopathies
Internal code and data layout for the project **shared_opposing_mechanisms_cmp**.

This repository documents the analysis pipeline and scripts used in the manuscript
"Leveraging the shared and opposing genetic mechanisms in the heritable cardiomyopathies".

All figures are generated using R (v4.3.1) and the tidyverse / Bioconductor ecosystem.  
A reproducible description of the R environment is provided in the **Code Availability** section.

Detailed documentation for specific figure inputs can be found in:

**`README_figures.md`**  

## Directory layout (high level)

- `code/` – all analysis scripts (R, Bash, etc.)
- `data_raw/` – raw input data (not for public sharing)
  - `data_raw/sumstats/` – original GWAS summary statistics (DCM, HCM)
- `data_processed/` – cleaned / harmonised data (to be defined later)
- `results/` – analysis outputs (to be defined later)

> Note: Only `code/` and small helper files will be shared publicly with the paper.
> Raw individual-level data and full GWAS summary statistics are subject to the
> original studies' data access policies.

---

## Step 1 – Obtain DCM and HCM GWAS summary statistics

1. Identify the two primary GWAS datasets:
   - DCM GWAS summary statistics from: https://www.nature.com/articles/s41588-024-01975-5
    Summary statistics for our GWAS meta-analyses have been made available for download through the Cardiovascular Disease Knowledge Portal (https://cvd.hugeamp.org/downloads.html)

   - HCM GWAS summary statistics from: https://www.nature.com/articles/s41588-025-02087-4
     Full GWAS summary statistics of HCM, HCMSARC−, HCMSARC+, MTAG and ten LV traits are available on the GWAS catalog (accession IDs GCST90435254 (https://www.ebi.ac.uk/gwas/studies/GCST90435254) –GCST90435267(https://www.ebi.ac.uk/gwas/studies/GCST90435267))


## Step 2 — Summary of QC for processed summary statistics (DCM/HCM)

QC filters:

1. Effect allele frequency filter:  
   - `EAFREQ >= 0.005 & EAFREQ <= 0.995`

2. Sample-size filters:  
   - For DCM GWAS: `N_cases >= 0.7 * max(N_cases)` 
   - For DCM MTAG: `N_Neff >= 0.7 * max(N_Neff)`
   - For HCM GWAS: `N_cases >= 0.96 * max(N_cases)`  
   - For HCM MTAG: `N_Neff >= 0.96 * max(N_Neff)` 

3. Exclusion of MHC-like extended region on chromosome 11:  
   - DCM: exclude variants with `CHR == 11 & POS ∈ [29,978,453; 80,288,956]`  
   - HCM: exclude variants with `CHR == 11 & POS ∈ [30,000,000; 80,000,000]`

Output: 
- harmonized_dcm.tsv.gz
- harmonized_hcm.tsv.gz
<!-- #### ===== packages and versions 
# conda activate LAVA_2024
# conda list
# conda env export > environment_LAVA_2024.yml
# R --version -->
<!-- Rscript -e "packageVersion('optparse')"
Rscript -e "packageVersion('data.table')"
Rscript -e "packageVersion('dplyr')"
Rscript -e "packageVersion('readr')"
Rscript -e "packageVersion('tidyr')" -->

---
## Step 3 - Genetic correlations

### 3.1 Global genetic correlation (rg)

Estimation of univariate SNP-heritability and bivariate genetic correlation using LD Score Regression (LDSC).

#### 3.1.1 LDSC software
We use LD Score Regression from Bulik-Sullivan et al.
Repository: https://github.com/bulik/ldsc

##### 3.1.2 Munging 

Example script:  
[`code/example_munging.sh`](code/example_munging.sh)

##### 3.1.3 LDSC Genetic Correlation 
Example script:  
[`code/example_ldsc_rg.sh`](code/example_ldsc_rg.sh)

##### TO ADD === version ==== source activate ldsc  # activate your conda env
Version recorded automatically in:
##### code/environment/ldsc_version.txt### 4.3 CC-GWAS: DCM vs HCM

##### 3.1.4 Figure 2a: Genetic correlation heatmap (DCM–HCM + MRI traits):
Example script:  
[`code/figure2a.r`](code/figure2a.R)


##### 3.1.4 Figure 2b: SNP effect concordance between DCM and HCM:
Example script:  
[`code/figure2b.r`](code/figure2a.R)

### 3.2 Local genetic correlation (LAVA)


LAVA was executed in a dedicated conda environment (`LAVA_2024`) to ensure consistent versions of PLINK, PLINK2, bcftools, and supporting R packages. Details and additional input files and script examples here: https://github.com/josefin-werme/LAVA 

Version used for our paper: https://github.com/josefin-werme/LAVA/releases/tag/v0.1.0 

#### 3.2.1 Activate environment
```bash
conda activate LAVA_2024
# R version (used for LAVA)
conda install -c bioconda r-base=4.4
```

```R
# Install LAVA from GitHub
install.packages("remotes")
remotes::install_github("josefin-werme/LAVA")
```

#### 3.2.2 Code example: [`code/run_lava_local_rg.r`](code/run_lava_local_rg.R)

#### 3.2.3 Figures 2d and 2f (LAVA locus annotation and Manhattan-type plots)


#### 3.3.2 Input Configuration File for LAVA

LAVA requires a simple tab-delimited configuration file specifying, for each phenotype: [`data/LAVA/input.info.txt`](data/input.info.txt)

## Step 4 – Installation of MTAG and CC-GWAS

MTAG and CC-GWAS are used for the case–case GWAS analyses in this project.

### 3.1 MTAG installation (Python 2.7)

MTAG is installed in a dedicated directory and conda environment:

- MTAG directory: `${HOME}/META_v2/MTAG`
- Conda environment: `mtag_py27` (Python 2.7)

The installation can be reproduced using:

```bash
bash code/environment/10_setup_mtag_ccgwas.sh
```

### 3.2 Inputs

This analysis requires harmonised DCM and HCM GWAS summary statistics as produced in Step 2 of the pipeline. All files must have consistent:

- chromosome (CHR)
- position (BP)
- SNP ID (SNP or rsID)
- effect allele (EA)
- non-effect allele (NEA)
- effect size (BETA)
- standard error (SE)
- p-value (P)
- effective sample size (Neff_p0.5 or N_cases depending on harmonisation)
    - `sum_stats/harmonized_dcm.tsv.gz`
    - `sum_stats/harmonized_hcm.tsv.gz`

Cardiac MRI traits for MTAG

These come from the cardiac MRI GWAS (Tadros et al.). Each file is formatted into MTAG-ready structure with:

- snpid
- a1, a2
- beta, se, pval
- freq
- n
    - `sum_stats/*_for_MTAG.txt` (Ecc, LVESVi, LVconc, etc.)

### === Sean processed? ====
#### 3.1 Run CC-GWAS

CC-GWAS contrasts the genetic architectures of DCM and HCM directly, modelling them as two case-control traits. This identifies variants with differential genetic effects.

Inputs for CC-GWAS:

- Harmonized DCM GWAS (with effective sample size Neff_p0.5)
- Harmonized HCM GWAS
- Prevalence estimates (population level + plausible bounds)
- LDSC-derived parameters:
    - SNP-heritability (DCM: ~0.142; HCM: ~0.1798)
    - DCM–HCM genetic correlation (≈ –0.5642)
    - Intercept (≈ 0.0107)
    - Number of causal SNPs (m ≈ 1200, based on typical polygenicity values)

The CC-GWAS function takes all these as explicit arguments.

```r
library(CCGWAS)

CCGWAS(
  outcome_file = "CC_GWAS__DCM__HCM.out",

  # Phenotype labels
  A_name = "DCM",
  B_name = "HCM",

  # Harmonised input sumstats
  sumstats_fileA1A0 = "sum_stats/Jurgens_DCM_LDSCinput_Neff.tsv.gz",
  sumstats_fileB1B0 = "sum_stats/Tadros_HCM_LDSCinput_Neff.tsv.gz",

  # Population prevalences (plausible range)
  K_A1A0      = 0.004,
  K_A1A0_high = 0.010,
  K_A1A0_low  = 0.002,

  K_B1B0      = 0.002,
  K_B1B0_high = 0.005,
  K_B1B0_low  = 0.001,

  # LDSC-derived genomic parameters
  h2l_A1A0           = 0.142,     # liability-scale h2 for DCM
  h2l_B1B0           = 0.1798,    # liability-scale h2 for HCM
  rg_A1A0_B1B0       = -0.5642,   # DCM–HCM genetic correlation
  intercept_A1A0_B1B0 = 0.0107,   # LDSC intercept (shared inflation)
  m                   = 1200,     # approximate causal SNP count

  # Effective sample sizes
  N_A1 = 9365  * 0.90,   # DCM cases (post-QC)
  N_A0 = 946368 * 0.90,  # DCM controls (post-QC)
  N_B1 = 5900   * 0.85,  # HCM cases (post-QC)
  N_B0 = 68359  * 0.85,  # HCM controls (post-QC)

  # Control overlap
  N_overlap_A0B0 = 40283 * 0.85
)
```
#### 3.2 Excluding MYBPC3 region (core snippet)

The MYBPC3 locus (chr11:30–80 Mb) is removed to avoid dominating downstream analyses due to extreme effect sizes in HCM.


#### 3.3 Preparing CC-GWAS for MTAG (description)

This step reformats the CC-GWAS output so it can be used as one of the traits in MTAG.

We:

- Keep core columns (snpid, a1, a2, beta, se, pval, z)
- Merge allele frequencies from DCM (EAFREQ)
- Compute minor allele frequency (MAF)
- Compute effective sample size per SNP (n = 1 / [2 * MAF * (1−MAF) * SE²])
- Export to MTAG format

Used later together with MRI traits in CC-MTAG.

## 4 Shared-effect meta analysis

2.5 Shared-effects meta-analysis of DCM and HCM

This step performs a shared-effects meta-analysis of DCM and HCM GWAS, under the assumption that the two traits share the same genetic effect (genetic correlation constrained to 1), while accounting for sample overlap. It uses MTAG with --equal_h2 and --perfect_gencov and then refines effect estimates with a secondary random-effects meta-analysis for loci with strong evidence.

2.5.1 Inputs

Two GWAS summary statistic files, pre-formatted for MTAG:
sum_stats/DCM_for_ccMTAG.txt
sum_stats/HCM_for_ccMTAG.txt

Each file contains (at least):
- snpid – SNP identifier (rsID or equivalent)
- A1, A2 – effect and non-effect allele
- beta, se, pval – effect size, SE, p-value
- freq – effect allele frequency
- n – effective sample size

These files are derived from the harmonised DCM and HCM GWAS, after standard QC and alignment.
## 4 Shared-effect meta analysis
## 4 Shared-effect meta analysis
## 4 Shared-effect meta analysis

## 5 Shared-effect meta analysis
