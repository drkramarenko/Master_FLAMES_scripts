#!/bin/bash
#SBATCH --job-name=ldsc_rg_bmi
#SBATCH --time=2:00:00
#SBATCH --mem=16G
source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh

conda activate ldsc  # Activate your conda environment if needed

WD_TOOLS=/home/dkramarenk/projects/tools
WD_PROJECT=/home/dkramarenk/projects/LAVA/DCM_HCM
INPUT_DIR=${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE
BMI_FILE=${INPUT_DIR}/bmi_clean.sumstats.gz

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

mkdir -p results_rg
mkdir -p logs

for BASE_NAME in "${FILES[@]}"; do
  echo "🔍 Processing: $BASE_NAME"
  BASE_NOEXT=${BASE_NAME%.sumstats.gz}
  OUT_DIR="results_rg/${BASE_NOEXT}"
  mkdir -p "$OUT_DIR"

  python ${WD_TOOLS}/ldsc/ldsc.py \
    --rg ${INPUT_DIR}/${BASE_NAME},${BMI_FILE} \
    --ref-ld-chr /gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/eur_w_ld_chr/ \
    --w-ld-chr /gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/eur_w_ld_chr/ \
    --out ${OUT_DIR}/ref_eur_bmi_meta
done