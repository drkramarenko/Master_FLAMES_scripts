# ================================
# LDSC master script
# ================================

# -------------------------------
# Creating an environment
# -------------------------------

# https://kevinlkx.github.io/analysis_pipelines/sldsc_pipeline.html

# cd /home/dkramarenk/projects/tools
# cd ldsc
# conda install -c bioconda htslib
# conda install -c bioconda csvtk
# conda env create --file environment.yml

# /home/dkramarenk/projects/tools/baselineLF_v2.2.UKB.tar.gz

conda activate ldsc

export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"
export WD_TOOLS="/home/dkramarenk/projects/tools/"

cd ${WD_PROJECT}
cd ${WD_TOOLS}/LDSC_files

srun -n 1 -c 16 -t 2:00:00 --pty /bin/bash 
# ================================
# Downloading other traits
# ================================

# ================================
# Munging HF_Henry_2025
# ================================
# HF_Henry_2025 https://www.nature.com/articles/s41588-024-02064-3#data-availability
cd /${WD_PROJECT}/data/sumst_processed/downloaded/HF_Henry_2025

INPUT_DIR=/${WD_PROJECT}/data/sumst_processed/downloaded/HF_Henry_2025
OUTPUT_DIR=/${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE
MERGE_ALLELES=/gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/w_hm3.snplist

for file in ${INPUT_DIR}/*ALL.tsv.gz; do
  base=$(basename "$file" .tsv.gz)
  echo "Processing: $base"

  ${WD_TOOLS}/ldsc/munge_sumstats.py \
    --out ${OUTPUT_DIR}/${base} \
    --merge-alleles ${MERGE_ALLELES} \
    --N-col N_total \
    --chunksize 500000 \
    --a1 A1 \
    --a2 A2 \
    --snp rsID \
    --sumstats "$file" \
    --p pval \
    --signed-sumstats A1_beta,0
done
# ================================
# Munging AF_Roselli_2025
# ================================
# AF_Roselli_2025 https://cvd.hugeamp.org/downloads.html#summary
cd /${WD_PROJECT}/data/sumst_processed/downloaded/AF_Roselli_2025

INPUT_DIR=/${WD_PROJECT}/data/sumst_processed/downloaded/AF_Roselli_2025
OUTPUT_DIR=/${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE
MERGE_ALLELES=/gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/w_hm3.snplist

zcat AF_GWAS_AFGenPlus_commonFreq_ALLv21.txt.gz | \
awk 'BEGIN {OFS="\t"}
NR==1 {
  for (i = 1; i <= NF; i++) {
    if ($i == "log(P)") logp_col = i
  }
  print $0, "p_val"
  next
}
{
  logp = $logp_col
  p = 10 ^ (logp)  # since log(P) = log10(P)
  print $0, p
}' > AF_GWAS_AFGenPlus_commonFreq_ALLv21_withP.txt

cut -f2,5- AF_GWAS_AFGenPlus_commonFreq_ALLv21_withP.txt > AF_for_munge.txt

for file in ${INPUT_DIR}/AF_for_munge.txt; do
  base=$(basename "$file" .tsv.gz)
  echo "Processing: $base"

  ${WD_TOOLS}/ldsc/munge_sumstats.py \
    --out ${OUTPUT_DIR}/${base} \
    --merge-alleles ${MERGE_ALLELES} \
    --N-col n_total \
    --chunksize 500000 \
    --a1 Allele1 \
    --a2 Allele2 \
    --snp rsid \
    --sumstats "$file" \
    --p p_val \
    --signed-sumstats Effect,0
done

# ================================
# Munging CAD_Aragam_2022
# ================================
cd /${WD_PROJECT}/data/sumst_processed/downloaded/CAD_Aragam_2022


zcat CAD_Aragam_2022.h.tsv.gz \
  | cut -f24,3,4,5,6,8,21,26 \
  | gzip > CAD_Aragam_2022_clean.tsv.gz

INPUT_DIR=/${WD_PROJECT}/data/sumst_processed/downloaded/CAD_Aragam_2022
OUTPUT_DIR=/${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE
MERGE_ALLELES=/gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/w_hm3.snplist

for file in ${INPUT_DIR}/CAD_Aragam_2022_clean.tsv.gz; do
  base=$(basename "$file" .tsv.gz)
  echo "Processing: $base"

  ${WD_TOOLS}/ldsc/munge_sumstats.py \
    --out ${OUTPUT_DIR}/${base} \
    --merge-alleles ${MERGE_ALLELES} \
    --N-col n \
    --chunksize 500000 \
    --a1 effect_allele \
    --a2 other_allele \
    --snp variant_id \
    --sumstats "$file" \
    --p p_value \
    --signed-sumstats beta,0
done

# ================================
# 1. S-LDSC: Partial heritability
# ================================
### 
# Convert GWAS summary stats to the LDSC .sumstats format using munge_sumstats.py
# ldsc wiki “Summary-Statistics-File-Format”
# Note: you may need to add an option --chunksize 500000 to munge_sumstats.py command

## Run ldsc on your GWAS summary statistics using baseline-LD model annotations and your new annotation
# Compute the annotation conditional on baselineLD model: controlling for the annotation categories of the 
# full baselineLD model, using a comma-separated list of annotation file prefixes with --ref-ld-chr

cd ${WD_PROJECT}//data/sumst_processed/ldsc/PART

# tutorial https://github.com/bulik/ldsc/wiki/Partitioned-Heritability
# https://kevinlkx.github.io/analysis_pipelines/sldsc_pipeline.html#prepare_gwas_summary_stats_in_ldsc_sumstats_format

# -------------------------------
# 1.1 S-LDSC: Partial heritability: Annotation file
# -------------------------------

# 1. Extract genome-wide significant SNPs from your sumstats
for name in CC_GWAS CC_MTAG; do
  awk '
  BEGIN {
    FS = OFS = "\t"
  }
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      if ($i == "CHR") chr = i
      if ($i == "BP") bp = i
      if ($i == "P") pval = i
      if ($i == "rsID") rsid = i
    }
    next
  }
  $pval < 5e-8 {
    start = $bp - 250000
    end = $bp + 250000
    if (start < 0) start = 0
    print "chr"$chr, start, end, $rsid
  }
  ' "/${WD_PROJECT}/data/sumst_processed/${name}_37_exclMYBPC3reg_cleaned.txt" > "gws_loci_${name,,}.bed"
done

cd ${WD_PROJECT}//data/sumst_processed/ldsc/PART

for name in CC_GWAS CC_MTAG; do
# /gpfs/work5/0/gusr0607/dkramarenko/tools/bedtools2/bin/bedtools = bedtools
/gpfs/work5/0/gusr0607/dkramarenko/tools/bedtools2/bin/bedtools sort -i "gws_loci_${name,,}.bed" | /gpfs/work5/0/gusr0607/dkramarenko/tools/bedtools2/bin/bedtools merge -i - > "gws_loci_${name,,}_merged.bed"
done

cat gws_loci_cc_gwas.bed gws_loci_cc_mtag.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  > gws_loci_cc_combined_merged.bed

cd ${WD_PROJECT}//data/sumst_processed/ldsc/PART

# 2. Generate .annot.gz

poloioio files using make_annot.py

for name in CC_GWAS CC_MTAG; do
for chr in {1..22}; do
Rscript ${WD_PROJECT}/scripts//make_ldsc_binary_annot.R \
  gws_loci_${name,,}_merged.bed \
  ${WD_TOOLS}/LDSC_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
  annot_shared_ld/gws_loci_${name,,}.${chr}.annot.gz "full-annot"
done
done
conda activate LAVA_2024 

for chr in {1..22}; do
Rscript ${WD_PROJECT}/scripts//make_ldsc_binary_annot.R \
  gws_loci_cc_combined_merged.bed \
  ${WD_TOOLS}/LDSC_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
  annot_shared_ld/gws_loci_cc_combined_merged.${chr}.annot.gz "full-annot"
done

for chr in {1..22}; do
python ${WD_TOOLS}/ldsc/ldsc.py \
  --l2 \
  --bfile ${WD_TOOLS}/LDSC_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
  --print-snps ${WD_TOOLS}/LDSC_files/hm3_no_MHC.list.txt \
  --ld-wind-cm 1 \
  --annot annot_shared_ld/gws_loci_cc_combined_merged.${chr}.annot.gz \
  --out annot_shared_ld/gws_loci_cc_combined_merged.${chr}
done

# 3. Making LD scores for annot file
for name in CC_GWAS CC_MTAG; do
for chr in {1..22}; do
python ${WD_TOOLS}/ldsc/ldsc.py \
  --l2 \
  --bfile ${WD_TOOLS}/LDSC_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
  --print-snps ${WD_TOOLS}/LDSC_files/hm3_no_MHC.list.txt \
  --ld-wind-cm 1 \
  --annot annot_shared_ld/gws_loci_${name,,}.${chr}.annot.gz \
  --out annot_shared_ld/gws_loci_${name,,}.${chr}
done
done 

# zcat  annot_shared_ld/gws_loci.1.l2.ldscore.gz | cut -f4 | sort | uniq -c
# zcat /${WD_TOOLS}/LDSC_files/baselineLD.1.l2.ldscore.gz | cut -f3 > baseline.snps
# zcat annot_shared/gws_loci.1.l2.ldscore.gz | cut -f3 > custom.snps

# diff baseline.snps custom.snps

# 4. Part. heritability for CC GWAS / CC MTAG

python ${WD_TOOLS}/ldsc/ldsc.py \
--h2 /${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE/DCM_GWAS_37_exclMYBPC3reg.sumstats.gz \
--ref-ld-chr /${WD_TOOLS}/LDSC_files/baselineLD.,/${WD_PROJECT}/data/sumst_processed/ldsc/PART/annot_shared_ld/gws_loci_cc_combined_merged. \
--frqfile-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot --print-cov --print-coefficients --print-delete-vals \
--out part_out_cc_gwas/cc_loci_DCM_sign_baselineLD

python ${WD_TOOLS}/ldsc/ldsc.py \
--h2 /${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE/HCM_GWAS_37_exclMYBPC3reg.sumstats.gz \
--ref-ld-chr /${WD_TOOLS}/LDSC_files/baselineLD.,/${WD_PROJECT}/data/sumst_processed/ldsc/PART/annot_shared_ld/gws_loci_cc_combined_merged. \
--frqfile-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot --print-cov --print-coefficients --print-delete-vals \
--out part_out_cc_gwas/cc_loci_HCM_sign_baselineLD


python ${WD_TOOLS}/ldsc/ldsc.py \
--h2 /${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE/shared_meta_analysis_DCM_HCM.sumstats.gz \
--ref-ld-chr /${WD_TOOLS}/LDSC_files/baselineLD.,/${WD_PROJECT}/data/sumst_processed/ldsc/PART/annot_shared_ld/gws_loci_cc_mtag. \
--frqfile-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot --print-cov --print-coefficients --print-delete-vals \
--out part_out_cc_mtag/cc_mtag_shared_sign_baselineLD

export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"
scp  dkramarenk@snellius.surf.nl:/${WD_PROJECT}//data/sumst_processed/ldsc/PART/part_out_cc_gwas/* /Users/drkramarenko/Library/CloudStorage/OneDrive-Personal/Computation/GWAS/Project_SHARED_DCM_HCM/ldsc


# ================================
# 2. Genetic correlation
# ================================

ll /gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/SUMSTATS/GWAS_sumstats /gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/SUMSTATS/GWAS_sumstats
export SJ_SUMST="/gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/MUNGE/TOTAL"

cd ${WD_PROJECT}//data/sumst_processed/ldsc/CORRELATION
# -------------------------------
# 3.1.1 Rg: all traits
# -------------------------------

sbatch ${WD_PROJECT}/scripts/eur_ref_ldsc_rg.sh
sbatch --export=BASE_NAME="shared_meta_analysis_DCM_HCM.sumstats.gz" ${WD_PROJECT}/scripts/eur_ref_ldsc_rg.sh
sbatch --export=BASE_NAME="DCM_GWAS_37_exclMYBPC3reg.sumstats.gz" ${WD_PROJECT}/scripts/eur_ref_ldsc_rg.sh
sbatch --export=BASE_NAME="DCM_MTAG_37_exclMYBPC3reg.sumstats.gz" ${WD_PROJECT}/scripts/eur_ref_ldsc_rg.sh
sbatch --export=BASE_NAME="HCM_GWAS_37_exclMYBPC3reg.sumstats.gz" ${WD_PROJECT}/scripts/eur_ref_ldsc_rg.sh
sbatch --export=BASE_NAME="HCM_MTAG_37_exclMYBPC3reg.sumstats.gz" ${WD_PROJECT}/scripts/eur_ref_ldsc_rg.sh
sbatch --export=BASE_NAME="CC_MTAG_37_exclMYBPC3reg.sumstats.gz" ${WD_PROJECT}/scripts/eur_ref_ldsc_rg.sh
sbatch --export=BASE_NAME="CC_GWAS_37_exclMYBPC3reg.sumstats.gz" ${WD_PROJECT}/scripts/eur_ref_ldsc_rg.sh
sbatch --export=BASE_NAME="FORMAT-METAL_Pheno1_ALL.sumstats.gz" ${WD_PROJECT}/scripts/eur_ref_ldsc_rg.sh
sbatch --export=BASE_NAME="FORMAT-METAL_Pheno2_ALL.sumstats.gz" ${WD_PROJECT}/scripts/eur_ref_ldsc_rg.sh
sbatch --export=BASE_NAME="FORMAT-METAL_Pheno3_ALL.sumstats.gz" ${WD_PROJECT}/scripts/eur_ref_ldsc_rg.sh
sbatch --export=BASE_NAME="FORMAT-METAL_Pheno4_ALL.sumstats.gz" ${WD_PROJECT}/scripts/eur_ref_ldsc_rg.sh
sbatch --export=BASE_NAME="AF_for_munge.txt.sumstats.gz" ${WD_PROJECT}/scripts/eur_ref_ldsc_rg.sh
sbatch --export=BASE_NAME="CAD_Aragam_2022_clean.sumstats.gz" ${WD_PROJECT}/scripts/eur_ref_ldsc_rg.sh

# -------------------------------
# 3.1.2 Rg:  + bmi
# -------------------------------

################################################## 
awk 'NR==1 || !/NA/' bmi.giant-ukbb.meta-analysis.combined.23May2018.txt > bmi_clean.txt

awk 'BEGIN {OFS="\t"} {gsub(/ +/, "\t"); print}' bmi_clean.txt > bmi_clean_t.txt

awk 'BEGIN {OFS="\t"} 
NR==1 {
  for (i = 1; i <= NF; i++) {
    if ($i == "SNP") snp_col = i
  }
  print
  next
}
{
  split($snp_col, parts, ":")
  $snp_col = parts[1]
  print
}' bmi_clean_t.txt > bmi_clean_rsID.txt

${WD_TOOLS}/ldsc/munge_sumstats.py \
--out /gpfs/work5/0/gusr0607/dkramarenko/LAVA/DCM_HCM/data/sumst_processed/ldsc/MUNGE/bmi_clean \
--merge-alleles /gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/w_hm3.snplist \
--N-col N \
--chunksize 500000 \
--a1 Tested_Allele \
--a2 Other_Allele \
--snp SNP \
--sumstats /gpfs/work5/0/gusr0607/dkramarenko/LAVA/DCM_HCM/data/sumst_processed/bmi_clean_rsID.txt \
--p P 

# python ${WD_TOOLS}/ldsc/ldsc.py \
#   --rg ${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE/shared_meta_analysis_DCM_HCM.sumstats.gz,/gpfs/work5/0/gusr0607/dkramarenko/LAVA/DCM_HCM/data/sumst_processed/ldsc/MUNGE/bmi_clean.sumstats.gz \
#   --ref-ld-chr ${WD_TOOLS}/LDSC_files/baselineLD. \
#   --w-ld-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
#   --out results_rg/shared_bmi_meta


# genetic_correlation_summary.R

# Set paths
WD_TOOLS=/home/dkramarenk/projects/tools
WD_PROJECT=/home/dkramarenk/projects/LAVA/DCM_HCM
INPUT_DIR=${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE
BMI_FILE=${INPUT_DIR}/bmi_clean.sumstats.gz

# List of files
FILES=(
  shared_meta_analysis_DCM_HCM.sumstats.gz
  DCM_GWAS_37_exclMYBPC3reg.sumstats.gz
  DCM_MTAG_37_exclMYBPC3reg.sumstats.gz
  HCM_GWAS_37_exclMYBPC3reg.sumstats.gz
  HCM_MTAG_37_exclMYBPC3reg.sumstats.gz
  CC_MTAG_37_exclMYBPC3reg.sumstats.gz
  CC_GWAS_37_exclMYBPC3reg.sumstats.gz
  FORMAT-METAL_Pheno1_ALL.sumstats.gz
  FORMAT-METAL_Pheno2_ALL.sumstats.gz
  FORMAT-METAL_Pheno3_ALL.sumstats.gz
  FORMAT-METAL_Pheno4_ALL.sumstats.gz
  AF_for_munge.txt.sumstats.gz
  CAD_Aragam_2022_clean.sumstats.gz
)

# Run loop
for BASE_NAME in "${FILES[@]}"; do
  echo "🔍 Processing: $BASE_NAME"

  BASE_NOEXT=${BASE_NAME%.sumstats.gz}
  OUT_DIR="results_rg/${BASE_NOEXT}"

  # Create output directory if it doesn't exist
  mkdir -p "$OUT_DIR"

  python ${WD_TOOLS}/ldsc/ldsc.py \
    --rg ${INPUT_DIR}/${BASE_NAME},${BMI_FILE} \
    --ref-ld-chr /gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/eur_w_ld_chr/ \
    --w-ld-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
    --out ${OUT_DIR}/bmi_meta
done

scp  dkramarenk@snellius.surf.nl:/${WD_PROJECT}//data/sumst_processed/ldsc/CORRELATION/results_rg/*.csv /Users/drkramarenko/Library/CloudStorage/OneDrive-Personal/Computation/GWAS/Project_SHARED_DCM_HCM/ldsc
scp  dkramarenk@snellius.surf.nl:/${WD_PROJECT}//data/sumst_processed/ldsc/CORRELATION/eur_ref_genetic_correlation_traits.csv /Users/drkramarenko/Library/CloudStorage/OneDrive-Personal/Computation/GWAS/Project_SHARED_DCM_HCM/ldsc
scp  dkramarenk@snellius.surf.nl:/${WD_PROJECT}//data/sumst_processed/ldsc/CORRELATION/results_rg_pairwise/eur_ref_results_rg_pairwise_summary.csv /Users/drkramarenko/Library/CloudStorage/OneDrive-Personal/Computation/GWAS/Project_SHARED_DCM_HCM/ldsc

# -------------------------------
# 3.1.2 Rg: pairwise
# -------------------------------

cd ${WD_PROJECT}//data/sumst_processed/ldsc/CORRELATION

# Activate LDSC environment
conda activate ldsc

# Set paths
WD_TOOLS=/home/dkramarenk/projects/tools
WD_PROJECT=/home/dkramarenk/projects/LAVA/DCM_HCM
INPUT_DIR=${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE
OUTPUT_DIR=${WD_PROJECT}/data/sumst_processed/ldsc/CORRELATION/results_rg_pairwise

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# List of summary statistics files
FILES=(
  shared_meta_analysis_DCM_HCM.sumstats.gz
  DCM_GWAS_37_exclMYBPC3reg.sumstats.gz
  DCM_MTAG_37_exclMYBPC3reg.sumstats.gz
  HCM_GWAS_37_exclMYBPC3reg.sumstats.gz
  HCM_MTAG_37_exclMYBPC3reg.sumstats.gz
  CC_MTAG_37_exclMYBPC3reg.sumstats.gz
  CC_GWAS_37_exclMYBPC3reg.sumstats.gz
  FORMAT-METAL_Pheno1_ALL.sumstats.gz
  FORMAT-METAL_Pheno2_ALL.sumstats.gz
  FORMAT-METAL_Pheno3_ALL.sumstats.gz
  FORMAT-METAL_Pheno4_ALL.sumstats.gz
  AF_for_munge.txt.sumstats.gz
  CAD_Aragam_2022_clean.sumstats.gz
)

# Loop over all unique pairs
for ((i=0; i<${#FILES[@]}; i++)); do
  for ((j=i+1; j<${#FILES[@]}; j++)); do
    FILE1=${FILES[$i]}
    FILE2=${FILES[$j]}

    BASE1=${FILE1%.sumstats.gz}
    BASE2=${FILE2%.sumstats.gz}
    
    echo "🔍 Calculating RG: $BASE1 vs $BASE2"

    python ${WD_TOOLS}/ldsc/ldsc.py \
      --rg ${INPUT_DIR}/${FILE1},${INPUT_DIR}/${FILE2} \
      --ref-ld-chr /gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/eur_w_ld_chr/ \
      --w-ld-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
      --out ${OUTPUT_DIR}/${BASE1}__vs__${BASE2}
  done
done

# ================================
# 2.2 Genetic correlation - rerun raw data
# ================================

### HCM
awk '
BEGIN {
  FS = OFS = "\t"
}
NR==1 {
  for (i = 1; i <= NF; i++) {
    col[$i] = i
  }
  print "rsid", "effect_allele", "noneffect_allele", "beta", "se", "pvalue", "n_samples"
  next
}
{
  print $col["rsid"], $col["effect_allele"], $col["noneffect_allele"], $col["beta"], $col["se"], $col["pvalue"], $col["n_samples"]
}
' hcm.gwama.250621.gel.nofin.txt > hcm_cleaned.txt


${WD_TOOLS}/ldsc/munge_sumstats.py \
--out /gpfs/work5/0/gusr0607/dkramarenko/LAVA/DCM_HCM/data/sumst_processed/ldsc/MUNGE/hcm_raw \
--merge-alleles /gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/w_hm3.snplist \
--N-col n_samples \
--chunksize 500000 \
--a1 effect_allele \
--a2 noneffect_allele \
--snp rsid \
--sumstats /gpfs/work5/0/gusr0607/dkramarenko/LAVA/DCM_HCM/data/raw/hcm_cleaned.txt \
--p pvalue 

### DCM
awk '
BEGIN {
  FS = OFS = "\t"
}
NR==1 {
  for (i = 1; i <= NF; i++) {
    col[$i] = i
  }
  print "rsID", "EA", "NEA", "BETA", "SE", "P", "EAFREQ", "N"
  next
}
{
  print $col["rsID"], $col["EA"], $col["NEA"], $col["BETA"], $col["SE"], $col["P"], $col["EAFREQ"], $col["N"]
}
' Garnier_Meder_Amsterdam_FinnGen_UKB_MGB__DCM__META1_chr1_22_MAF0.005.tsv > dcm_cleaned.txt

${WD_TOOLS}/ldsc/munge_sumstats.py \
--out /gpfs/work5/0/gusr0607/dkramarenko/LAVA/DCM_HCM/data/sumst_processed/ldsc/MUNGE/dcm_raw \
--merge-alleles /gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/w_hm3.snplist \
--N-col N \
--chunksize 500000 \
--a1 EA \
--a2 NEA \
--snp rsID \
--sumstats /gpfs/work5/0/gusr0607/dkramarenko/LAVA/DCM_HCM/data/raw/dcm_cleaned.txt \
--p pvalue 

## correlate
cd ${WD_PROJECT}//data/sumst_processed/ldsc/CORRELATION

WD_TOOLS=/home/dkramarenk/projects/tools
WD_PROJECT=/home/dkramarenk/projects/LAVA/DCM_HCM
INPUT_DIR=${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE
OUTPUT_DIR=${WD_PROJECT}/data/sumst_processed/ldsc/CORRELATION/results_rg_pairwise

FILE1=dcm_raw.sumstats.gz
FILE2=hcm_raw.sumstats.gz

BASE1=${FILE1%.sumstats.gz}
BASE2=${FILE2%.sumstats.gz}

echo "🔍 Calculating RG: $BASE1 vs $BASE2"

python ${WD_TOOLS}/ldsc/ldsc.py \
  --rg ${INPUT_DIR}/${FILE1},${INPUT_DIR}/${FILE2} \
  --ref-ld-chr /gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/eur_w_ld_chr/ \
  --w-ld-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
  --out ${OUTPUT_DIR}/${BASE1}__vs__${BASE2}

# ================================
# 3. Coloc
# ================================

# conda create -n coloc
# environment location: /home/dkramarenk/.conda/envs/coloc
# conda install -c bioconda r-base=4.3.2

export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"
export WD_TOOLS="/home/dkramarenk/projects/tools/"

cd ${WD_PROJECT}

conda activate coloc
# make a shared file
#------------
library(dplyr)
library(readr)

# Read data
df <- read_tsv("shared_meta_analysis_DCM_HCM.txt")

# Get top SNP by harmonized_pval
top_snp <- df %>%
  filter(!is.na(harmonized_pval)) %>%
  arrange(harmonized_pval) %>%
  slice(1)

chr <- top_snp$CHR
pos <- top_snp$BP
region <- df %>%
  filter(CHR == chr, BP >= (pos - 5e5), BP <= (pos + 5e5))

write_tsv(region, "shared_top_region_500kb.txt")
#------------
awk -F'\t' 'NR > 1 { print $5, $6 }' shared_top_region_500kb.txt > shared_top_coords.bed
awk '{start = $2 - 500000; if (start < 0) start = 0; end = $2 + 500000; print "chr"$1, start, end}' OFS='\t' shared_top_coords.bed > shared_top_region_500kb.bed

zcat GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Heart_Left_Ventricle.allpairs.txt.gz \
| awk 'BEGIN{OFS="\t"} NR>1 {
    split($2, a, "_");
    print a[1], a[2], $0
}' > gtex_chr_pos.tsv

awk '{print $1, $2, $2+1, $0}' OFS='\t' gtex_chr_pos.tsv > gtex_eqtl.bed

/gpfs/work5/0/gusr0607/dkramarenko/tools/bedtools2/bin/bedtools intersect -a gtex_eqtl.bed -b shared_top_region_500kb.bed > gtex_filtered.bed


# coloc
#------------
library(coloc)
library(dplyr)
library(readr)
library(data.table)
library(tidyr)
library(stringr)

gwas_df <- read_tsv("shared_top_region_500kb.txt") 
# Filter out NAs just in case
gwas_df1 <- gwas_df %>% 
filter(!is.na(harmonized_beta), !is.na(harmonized_se), !is.na(harmonized_pval)) %>%
  transmute(
        CHR, BP, A1, A2,
         beta = harmonized_beta,
         se = harmonized_se,
         pval = harmonized_pval,
         MAF = meta_freq,
         position = BP,
         chr = CHR) %>%
 mutate(snp = paste0("chr", CHR, "_", BP, "_", A2, "_", A1))

gwas_df2 <- gwas_df %>%
  transmute(
        CHR, BP, A1.1 = A2, A2.2 = A1,
         beta = harmonized_beta*-1,
         se = harmonized_se,
         pval = harmonized_pval,
         MAF = meta_freq,
         position = BP,
         chr = CHR) %>%
  rename(A1 = A1.1, A2  = A2.2) %>%
 mutate(snp = paste0("chr", CHR, "_", BP, "_", A2, "_", A1))

gwas_both <- bind_rows(gwas_df1,gwas_df2)
qtl_raw <- fread("gtex_filtered.bed") 
colnames(qtl_raw) <- c("chr", "pos", "pos2", "chr_temp", "pos_temp", "gene_id", "variant_id", "tss_distance", "ma_samples",
                    "ma_count", "maf", "pval_nominal", "slope", "slope_se")

qtl_raw$snp <- sub("_b38$", "", qtl_raw$variant_id)
  
common_snps <- intersect(gwas_both$snp, qtl_raw$snp)

gwas_sub <- gwas_both[gwas_both$snp %in% common_snps, ]
eqtl_sub <- qtl_raw[qtl_raw$snp %in% common_snps, ]

# Sort to align
gwas_sub <- gwas_sub[order(gwas_sub$snp), ]
eqtl_sub <- eqtl_sub[order(eqtl_sub$snp), ]

eqtl_sub <- eqtl_sub %>%
  group_by(variant_id) %>%
  slice(1) %>%
  ungroup()

dataset1 <- list(
  beta = gwas_sub$beta,
  varbeta = gwas_sub$se^2,
  snp = gwas_sub$snp,
  #MAF = gwas_sub$MAF,
  type = "cc",
  N = 27388  # adjust to your sample size
)

dataset2 <- list(
  beta = eqtl_sub$slope,
  varbeta = eqtl_sub$slope_se^2,
  snp = sub("_b38$", "", eqtl_sub$variant_id),
  #MAF = eqtl_sub$maf,
  type = "cc",
  N = 1000   # adjust to GTEx sample size for tissue
)

result <- coloc.abf(dataset1, dataset2)

PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
 3.33e-04  8.64e-01  3.88e-05  1.01e-01  3.51e-02 
[1] "PP abf for shared variant: 3.51%"

PP.H0 = 0.0003   # No association with either trait
PP.H1 = 86.4%    # Only the GWAS (cardiomyopathy) has a signal
PP.H2 = 0.004%   # Only the eQTL has a signal
PP.H3 = 10.1%    # Both have signals, but from *different* variants
PP.H4 = 3.5%     # Both share the *same* causal variant
In check_dataset(d = dataset2, 2) : minimum p value is: 0.0022944


# ================================
# 5. MRI traits
# ================================


cd /archive/expcard/Lisa_Projects/HCM-GWAS/HCM-2022/

ll /${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE

# Process MRI traits

# ---------------------------------------
#!/bin/bash
#SBATCH --job-name=LV_traits
#SBATCH --time=02:00:00
#SBATCH --mem=32G

source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh
conda activate ldsc

set -euo pipefail

# Define directories
INPUT_DIR="/home/dkramarenk/projects/LAVA/DCM_HCM/data/sumst_processed"
OUTPUT_DIR="/home/dkramarenk/projects/LAVA/DCM_HCM/data/sumst_processed"

# Create output directory if needed

# Loop over each Jan2022 .csv.gz file
for file in "$INPUT_DIR"/*Jan2022*.csv.gz; do
  base=$(basename "$file" .csv.gz)
  output_file="${OUTPUT_DIR}/${base}.tsv.gz"

  echo "🔄 Converting and cleaning: $file"

  # Convert CSV to TSV, remove rows with any NA, and save as .tsv.gz
  zcat "$file" | \
    awk -F',' 'NR==1 || ($0 !~ /NA/)' | \
    tr ',' '\t' | \
    gzip > "$output_file"

  echo "✅ Saved cleaned file: $output_file"
done

# ---------------------------------------

#!/bin/bash
#SBATCH --job-name=LV_traits
#SBATCH --time=02:00:00
#SBATCH --mem=32G

source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh
conda activate ldsc

WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"
WD_TOOLS="/home/dkramarenk/projects/tools"

INPUT_DIR="/home/dkramarenk/projects/LAVA/DCM_HCM/data/sumst_processed"
OUTPUT_DIR="${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE"
MERGE_ALLELES="/gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/w_hm3.snplist"

for file in ${INPUT_DIR}/*Jan2022*tsv.gz; do
# for file in ${INPUT_DIR}/maxWT_Jan2022_fullCADexcl.tsv.gz; do
  echo "Checking: $file"
  if [[ ! -f "$file" ]]; then
    echo "❌ Skipping missing file: $file"
    continue
  fi

  base=$(basename "$file" .csv.gz)
  echo "🔄 Processing: $base"

  python ${WD_TOOLS}/ldsc/munge_sumstats.py \
    --out "${OUTPUT_DIR}/${base}" \
    --merge-alleles "$MERGE_ALLELES" \
    --chunksize 500000 \
    --a1 ALLELE1 \
    --a2 ALLELE0 \
    --snp SNP \
    --p P_LINREG \
    --signed-sumstats BETA,0 \
    --frq A1FREQ \
    --N '36083' \
    --sumstats "$file"
done

#### genetic correlation
# ---------------------------------------

#!/bin/bash
#SBATCH --job-name=LV_traits
#SBATCH --time=02:00:00
#SBATCH --mem=32G

source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh
conda activate ldsc

WD_TOOLS=/home/dkramarenk/projects/tools
WD_PROJECT=/home/dkramarenk/projects/LAVA/DCM_HCM
INPUT_DIR=${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE
OUTPUT_DIR=${WD_PROJECT}/data/sumst_processed/ldsc/CORRELATION/results_rg_gwas

# List of files
FILES=(
  shared_meta_analysis_DCM_HCM.sumstats.gz
  DCM_GWAS_37_exclMYBPC3reg.sumstats.gz
  HCM_GWAS_37_exclMYBPC3reg.sumstats.gz
  CC_MTAG_37_exclMYBPC3reg.sumstats.gz
  CC_GWAS_37_exclMYBPC3reg.sumstats.gz
  Ecc_global_Jan2022_fullCADexcl.tsv.gz.sumstats.gz  
  LVEF_Jan2022_fullCADexcl.tsv.gz.sumstats.gz
  Ell_global_Jan2022_fullCADexcl.tsv.gz.sumstats.gz  
  LVESVi_Jan2022_fullCADexcl.tsv.gz.sumstats.gz
  Err_global_Jan2022_fullCADexcl.tsv.gz.sumstats.gz  
  LVMi_Jan2022_fullCADexcl.tsv.gz.sumstats.gz
  LVconc_Jan2022_fullCADexcl.tsv.gz.sumstats.gz      
  meanWT_Jan2022_fullCADexcl.tsv.gz.sumstats.gz
  LVEDVi_Jan2022_fullCADexcl.tsv.gz.sumstats.gz
)

for ((i=0; i<${#FILES[@]}; i++)); do
  for ((j=i+1; j<${#FILES[@]}; j++)); do
    FILE1=${FILES[$i]}
    FILE2=${FILES[$j]}

    BASE1=${FILE1%.sumstats.gz}
    BASE2=${FILE2%.sumstats.gz}
    
    echo "🔍 Calculating RG: $BASE1 vs $BASE2"

    python ${WD_TOOLS}/ldsc/ldsc.py \
      --rg ${INPUT_DIR}/${FILE1},${INPUT_DIR}/${FILE2} \
      --ref-ld-chr /gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/eur_w_ld_chr/ \
      --w-ld-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
      --out ${OUTPUT_DIR}/${BASE1}__vs__${BASE2}
  done
done

# ---------------------------------------

#!/bin/bash
#SBATCH --job-name=LV_traits
#SBATCH --time=01:00:00
#SBATCH --mem=32G

source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh
conda activate ldsc

WD_TOOLS=/home/dkramarenk/projects/tools
WD_PROJECT=/home/dkramarenk/projects/LAVA/DCM_HCM

INPUT_DIR=${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE
FIRST_FILE=${INPUT_DIR}/maxWT_Jan2022_fullCADexcl.tsv.gz.sumstats.gz
OUTPUT_DIR=${WD_PROJECT}/data/sumst_processed/ldsc/CORRELATION/results_rg_gwas

# List of files
FILES=(
  shared_meta_analysis_DCM_HCM.sumstats.gz
  DCM_GWAS_37_exclMYBPC3reg.sumstats.gz
  HCM_GWAS_37_exclMYBPC3reg.sumstats.gz
  CC_MTAG_37_exclMYBPC3reg.sumstats.gz
  CC_GWAS_37_exclMYBPC3reg.sumstats.gz
  Ecc_global_Jan2022_fullCADexcl.tsv.gz.sumstats.gz  
  LVEF_Jan2022_fullCADexcl.tsv.gz.sumstats.gz
  Ell_global_Jan2022_fullCADexcl.tsv.gz.sumstats.gz  
  LVESVi_Jan2022_fullCADexcl.tsv.gz.sumstats.gz
  Err_global_Jan2022_fullCADexcl.tsv.gz.sumstats.gz  
  LVMi_Jan2022_fullCADexcl.tsv.gz.sumstats.gz
  LVconc_Jan2022_fullCADexcl.tsv.gz.sumstats.gz      
  meanWT_Jan2022_fullCADexcl.tsv.gz.sumstats.gz
  LVEDVi_Jan2022_fullCADexcl.tsv.gz.sumstats.gz
)

# Run loop
for BASE_NAME in "${FILES[@]}"; do
  echo "🔍 Processing: $BASE_NAME"

  python ${WD_TOOLS}/ldsc/ldsc.py \
    --rg ${INPUT_DIR}/${BASE_NAME},${FIRST_FILE} \
    --ref-ld-chr /gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/eur_w_ld_chr/ \
    --w-ld-chr ${WD_TOOLS}/LDSC_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
    --out ${OUTPUT_DIR}/${BASE_NAME}__vs__${FIRST_FILE}

done









scp  dkramarenk@snellius.surf.nl:/${WD_PROJECT}//data/sumst_processed/ldsc/CORRELATION/results_rg_gwas/LV_traits_pairwise_summary.csv /Users/drkramarenko/Library/CloudStorage/OneDrive-Personal/Computation/GWAS/Project_SHARED_DCM_HCM/ldsc

scp  dkramarenk@snellius.surf.nl://gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/eur_w_ld_chr/*.gz  /Users/drkramarenko/Library/CloudStorage/OneDrive-Personal/Computation/LDSC
scp  dkramarenk@snellius.surf.nl://gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/eur_w_ld_chr/w_hm3.snplist  /Users/drkramarenko/Library/CloudStorage/OneDrive-Personal/Computation/LDSC/eur_w_ld_chr
scp  dkramarenk@snellius.surf.nl://gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/eur_w_ld_chr/README  /Users/drkramarenko/Library/CloudStorage/OneDrive-Personal/Computation/LDSC/eur_w_ld_chr


HEARTFAIL.gwas.imputed_v3.both_sexes.with_major.tsv.bgz

scp  dkramarenk@snellius.surf.nl:/home/dkramarenk/test/HEARTFAIL.gwas.imputed_v3.both_sexes.with_major.tsv.bgz /Users/drkramarenko/Library/CloudStorage/OneDrive-Personal/Computation/GWAS/test
