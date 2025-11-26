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
