#!/bin/bash
#SBATCH --job-name=eur_ref_ldsc_rg
#SBATCH --array=0-65
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
source /sw/arch/RHEL8/EB_production/2022/software/Anaconda3/2022.05/etc/profile.d/conda.sh

conda activate ldsc  # activate your conda env

# === Paths ===
WD_TOOLS=/home/dkramarenk/projects/tools
WD_PROJECT=/home/dkramarenk/projects/LAVA/DCM_HCM
SJ_SUMST=/gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/MUNGE/TOTAL  

# === Dataset file to test against ===
# BASE_NAME="shared_meta_analysis_DCM_HCM.sumstats.gz"
BASE_PATH="${WD_PROJECT}/data/sumst_processed/ldsc/MUNGE/${BASE_NAME}"

# === Get Trait from Array Index ===
trait_file=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ${SJ_SUMST}/traitsall.txt)
trait=$(basename "$trait_file" .tsv.gz.sumstats.gz)
trait=${trait#GWAS_sumstats_EUR__invnorm_}
trait=${trait%__TOTALsample}

echo "Running LDSC RG for trait: $trait"

# Define output directory
OUT_DIR="results_rg/${BASE_NAME%.sumstats.gz}"

# Create it if it doesn't exist
if [ ! -d "$OUT_DIR" ]; then
  mkdir -p "$OUT_DIR"
fi

# Run LDSC
python ${WD_TOOLS}/ldsc/ldsc.py \
  --rg ${BASE_PATH},${SJ_SUMST}/${trait_file} \
  --ref-ld-chr $/gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/eur_w_ld_chr/ \
  --w-ld-chr /gpfs/work5/0/gusr0607/dominicz/GWAS_SCA/HeritabilityCorrelation/ldsc/eur_w_ld_chr/ \
  --out ${OUT_DIR}/eur_ref_${BASE_NAME%.sumstats.gz}_${trait}
