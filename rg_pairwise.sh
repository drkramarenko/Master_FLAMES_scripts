#!/bin/bash
#SBATCH --job-name=rg_pairwise
#SBATCH --time=2:00:00
#SBATCH --mem=16G

source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh

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
  HCM_GWAS_37_exclMYBPC3reg.sumstats.gz
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
      --w-ld-chr /gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/eur_w_ld_chr/ \
      --out ${OUTPUT_DIR}/eur_ref_${BASE1}__vs__${BASE2}
  done
done
