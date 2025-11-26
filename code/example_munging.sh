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