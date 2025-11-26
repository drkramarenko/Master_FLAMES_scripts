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

Code examples:
##### 3.1.2 Munging 
``` bash
# Define paths
PROJECT_DIR=/path/to/project                                        # Main project directory
LDSC_DIR=/path/to/ldsc                                              # LDSC installation directory
INPUT_FILE=${PROJECT_DIR}/data/sumstats_raw/<FILENAME>.tsv.gz       # Input file (example; replace <FILENAME> as needed)
OUTPUT_PREFIX=${PROJECT_DIR}/data/sumstats_ldsc_ready/<FILENAME>    # Output prefix (example; replace <FILENAME> as needed)
MERGE_ALLELES=${LDSC_DIR}/w_hm3.snplist                             # Merge-alleles HapMap3 file (for LDSC)

# Run LDSC munging
${LDSC_DIR}/munge_sumstats.py \
    --sumstats "${INPUT_FILE}" \
    --out "${OUTPUT_PREFIX}" \
    --merge-alleles "${MERGE_ALLELES}" \
    --snp rsID \
    --a1 A1 \
    --a2 A2 \
    --N-col N_total \
    --p pval \
    --signed-sumstats A1_beta,0 \
    --chunksize 500000
```

##### 3.1.3 LDSC Genetic Correlation 
``` bash
# Define paths
TOOLS_DIR=/path/to/tools                     # directory with ldsc/ and reference files
PROJECT_DIR=/path/to/project                 # project root directory
SUMST_DIR=/path/to/sumstats_ldsc_ready       # directory with munged sumstats
TRAIT_FILE=/path/to/other/trait.sumstats.gz  # second trait for comparison

BASE_NAME=<BASE_TRAIT>.sumstats.gz
BASE_PATH=${SUMST_DIR}/${BASE_NAME}

# 
OUT_DIR=${PROJECT_DIR}/results/ldsc_rg/<BASE_TRAIT>/
mkdir -p "${OUT_DIR}"

# === Run LDSC Genetic Correlation ===
python ${TOOLS_DIR}/ldsc/ldsc.py \
    --rg ${BASE_PATH},${TRAIT_FILE} \
    --ref-ld-chr ${TOOLS_DIR}/ldsc_ref/baselineLD. \
    --w-ld-chr ${TOOLS_DIR}/ldsc_ref/weights_hm3_noMHC. \
    --out ${OUT_DIR}/${BASE_NAME%.sumstats.gz}_vs_<TRAIT_NAME>
```

##### === version ==== source activate ldsc  # activate your conda env
Version recorded automatically in:
##### code/environment/ldsc_version.txt### 4.3 CC-GWAS: DCM vs HCM

##### 3.1.4 Figure 2a: Genetic correlation heatmap (DCM–HCM + MRI traits):
```R
# Input: ST3 Genetic correlations: DCM, HCM & cardiac MRI traits
# Output: Figure 2a

# Package versions used:
# readxl   1.4.5
# dplyr    1.1.4
# reshape2 1.4.4
# ggplot2  3.5.2

library(readxl)    # read_excel()
library(dplyr)     # %>%, mutate(), across()
library(reshape2)  # melt()
library(ggplot2)   # ggplot2 plotting

# --- Paths and input file ---
PROJECT_DIR <- "/path/to/project"
FIG_INPUT_DIR <- file.path(PROJECT_DIR, "figures/for_figures")
FIG_OUTPUT_DIR <- file.path(PROJECT_DIR, "figures/output")

HEATMAP_XLSX <- file.path(FIG_INPUT_DIR, "heatmap_main.xlsx")
HEATMAP_SHEET <- "DCM_HCM"

# --- Read and prepare data ---
gen_cor_DCM_HCM <- read_excel(HEATMAP_XLSX, sheet = HEATMAP_SHEET) %>%
  mutate(across(-rows, as.numeric))

# Long format for ggplot
long_data <- melt(
  gen_cor_DCM_HCM,
  id.vars = "rows",
  variable.name = "Variable",
  value.name = "Correlation"
)

# Ensure rows and columns use the same order
order_levels <- c("DCM", "Ecc", "LVESVi", "LVconc", "HCM")
long_data$rows     <- factor(long_data$rows,     levels = order_levels)
long_data$Variable <- factor(long_data$Variable, levels = order_levels)

# --- Create heatmap ---
DCM_HCM_mat <- ggplot(long_data, aes(x = Variable, y = rows, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 4) +
  scale_fill_gradient2(
    low  = "#078eff",   # blue
    mid  = "#f2f7f3",   # neutral
    high = "#B71C1C",   # red
    midpoint = 0,
    limits   = c(-1, 1)
  ) +
  labs(title = "", x = "", y = "") +
  theme_minimal()

print(DCM_HCM_mat)

# --- Save outputs ---
dir.create(FIG_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(FIG_OUTPUT_DIR, "DCM_HCM_mat.png"),
  plot     = DCM_HCM_mat,
  width    = 6,
  height   = 4,
  dpi      = 450,
  bg       = "transparent"
)

ggsave(
  filename = file.path(FIG_OUTPUT_DIR, "DCM_HCM_mat.pdf"),
  plot     = DCM_HCM_mat,
  width    = 6,
  height   = 4,
  device   = "pdf",
  bg       = "transparent"
)
```

##### 3.1.4 Figure 2b: SNP effect concordance between DCM and HCM:

```R
# Input: ST4 Comparison of lead variant effects: DCM vs HCM
# Output: Figure 2b
# ggplot2 3.5.2
# dplyr   1.1.4   (if used for preprocessing)
# data.table / readr – optional depending on how loci_for_mr_selected is loaded

library(ggplot2)

# --- Input object ---
# loci_for_mr_selected must contain:
#   DCM_BETA, HCM_BETA, DCM_SE, HCM_SE, SNP_label
# Load your data here if needed:
# loci_for_mr_selected <- readRDS(file.path(PROJECT_DIR, "data/processed/loci_for_mr_selected.rds"))

# --- Plot ---
effect_concordance_plot <- ggplot(
  loci_for_mr_selected,
  aes(x = DCM_BETA, y = HCM_BETA, color = SNP_label)
) +
  # Regression line
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "gray40",
    linetype = "dashed",
    linewidth = 0.5,
    fullrange = TRUE
  ) +
  # Points
  geom_point(alpha = 0.7, size = 3) +
  
  # Manual color scale
  scale_color_manual(
    values = c(
      "DCM"  = "#7030A0",
      "HCM"  = "#00B050",
      "Both" = "red4"
    )
  ) +
  
  # Error bars
  geom_errorbar(
    aes(ymin = HCM_BETA - HCM_SE, ymax = HCM_BETA + HCM_SE),
    width = 0,
    alpha  = 0.5
  ) +
  geom_errorbarh(
    aes(xmin = DCM_BETA - DCM_SE, xmax = DCM_BETA + DCM_SE),
    height = 0,
    alpha  = 0.5
  ) +

  # Reference lines
  geom_vline(xintercept = 0, color = "gray60") +
  geom_hline(yintercept = 0, color = "gray60") +
  
  # Labels
  labs(
    x = "DCM Beta",
    y = "HCM Beta",
    title = "",
    color = "Locus label"
  ) +
  
  theme_classic() +
  xlim(-0.4, 0.4) +
  ylim(-0.4, 0.55) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title  = element_text(size = 10),
    axis.text   = element_text(size = 7)
  )

# Display plot
print(effect_concordance_plot)

# --- Save outputs ---
ggsave(
  filename = file.path(FIG_OUTPUT_DIR, "figure2b_effect_concordance.png"),
  plot     = effect_concordance_plot,
  width    = 6,
  height   = 4,
  dpi      = 450,
  bg       = "transparent"
)

ggsave(
  filename = file.path(FIG_OUTPUT_DIR, "figure2b_effect_concordance.pdf"),
  plot     = effect_concordance_plot,
  width    = 6,
  height   = 4,
  device   = "pdf",
  bg       = "transparent"
)

```


### 3.2 Local genetic correlation (LAVA)


LAVA was executed in a dedicated conda environment (`LAVA_2024`) to ensure consistent versions of PLINK, PLINK2, bcftools, and supporting R packages.

#### Activate environment
```bash
conda activate LAVA_2024
# PLINK 1.9
conda install -c bioconda plink=1.90b6.21

# PLINK 2.0 (build compatible with LAVA)
conda install -c bioconda plink2=2.00a2.3
# (Note: newer PLINK2 builds exist but conflict with R requirements.)

# bcftools
conda install -c bioconda bcftools=1.17

# zstd (compression library sometimes required by PLINK2)
conda install -c bioconda zstd=1.5.5

# R version (used for LAVA)
conda install -c bioconda r-base=4.4
```
```R
# Install BiocManager if needed
install.packages("BiocManager")


# Required for genotype matrix handling
BiocManager::install("snpStats")

# Install LAVA from GitHub
install.packages("remotes")
remotes::install_github("josefin-werme/LAVA")
```

#### 3.3.2 Input Configuration File for LAVA

LAVA requires a simple tab-delimited configuration file specifying, for each phenotype:

1. phenotype name (string)
2. number of cases
3. number of controls
4. path to processed GWAS summary statistics (LAVA-formatted)

This file is used by the LAVA wrapper scripts to:
- detect which phenotypes to include,
- locate the pre-processed input summary statistics,
- determine sample sizes for univariate and bivariate tests.

### Example: `data/input.info.txt`
```text
phenotype   cases   controls    filename
DCM         37083   37083       /home/dkramarenk/projects/LAVA/DCM_HCM/data/sumst_processed/DCM_GWAS_37_exclMYBPC3reg.txt
HCM         21708   21708       /home/dkramarenk/projects/LAVA/DCM_HCM/data/sumst_processed/HCM_GWAS_37_exclMYBPC3reg.txt
```

---

### `scripts/run_lava_local_rg.R` (clean template)

```r
#!/usr/bin/env Rscript

## LAVA local rg analysis (univariate + bivariate) for a block of loci
## Usage:
## Rscript run_lava_local_rg.R \
##   <ref_prefix> <locfile> <input.info> <sample.overlap> <phenos> <out_prefix> <start> <stop>
##
## Example:
## Rscript run_lava_local_rg.R \
##   data/refdat/g1000_eur.maf01.ids \
##   data/blocks_s2500_m25_f1_w200.locfile \
##   data/input.info.txt \
##   data/sample.overlap.txt \
##   DCM;HCM \
##   results/lava_dcm_hcm_out \
##   1 5

suppressPackageStartupMessages({
  library(parallel)
  library(LAVA)
  library(snpStats)
  library(data.table)
  library(matrixsampling)
  library(keep)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 8) {
  stop(
    "Expected 8 arguments:\n",
    " 1) ref_prefix (e.g. data/refdat/g1000_eur.maf01.ids)\n",
    " 2) locfile (e.g. data/blocks_s2500_m25_f1_w200.locfile)\n",
    " 3) input.info (e.g. data/input.info.txt)\n",
    " 4) sample.overlap (e.g. data/sample.overlap.txt)\n",
    " 5) phenos (e.g. DCM;HCM)\n",
    " 6) out_prefix (e.g. results/lava_dcm_hcm_out)\n",
    " 7) start_locus (integer)\n",
    " 8) stop_locus (integer)\n",
    call. = FALSE
  )
}

ref_prefix   <- args[1]
locfile      <- args[2]
input_info   <- args[3]
overlap_file <- args[4]
phenos       <- unlist(strsplit(args[5], ";"))
out_prefix   <- args[6]
start_locus  <- as.integer(args[7])
stop_locus   <- as.integer(args[8])

if (any(is.na(c(start_locus, stop_locus)))) {
  stop("start_locus and stop_locus must be integers.", call. = FALSE)
}

univ_thresh <- 0.027  # example: 0.05 / N_loci (adjust as needed)

## ---- Read loci definition ----
loci <- LAVA::read.loci(locfile)
n_loci_total <- nrow(loci)

if (stop_locus > n_loci_total) {
  stop("stop_locus exceeds number of loci in locfile: ", n_loci_total, call. = FALSE)
}

message("Running LAVA for loci ", start_locus, " to ", stop_locus, " out of ", n_loci_total)

## ---- Process input (GWAS + LD + sample overlap) ----
input <- LAVA::process.input(
  input.info.file    = input_info,
  sample.overlap.file = overlap_file,
  ref.prefix          = ref_prefix,
  phenos              = phenos
)

## ---- Analysis loop ----
loc_indices <- start_locus:stop_locus
progress <- ceiling(stats::quantile(loc_indices, seq(0.05, 1, 0.05)))

univ_results  <- NULL
bivar_results <- NULL

message(Sys.time(), "  Starting LAVA analysis for ", length(loc_indices), " loci")

for (i in loc_indices) {

  if (i %in% progress) {
```

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
