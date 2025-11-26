#!/usr/bin/env bash
# 10_setup_mtag_ccgwas.sh
#
# Install MTAG (Python 2.7) and the CC-GWAS R package.
# This script is meant as REPRODUCIBLE DOCUMENTATION of the steps;
# depending on your cluster/conda setup, you may run commands
# interactively or adapt paths.
#
# MTAG will be installed in:  ${HOME}/META_v2/MTAG
# Conda env name:             mtag_py27
#
# Usage:
#   bash code/environment/10_setup_mtag_ccgwas.sh
#
# NOTE: Python 2.7 is end-of-life; this environment is isolated and
#       used only for MTAG.

set -euo pipefail

# ---------- Configuration ----------
MTAG_ROOT="${HOME}/META_v2"
MTAG_DIR="${MTAG_ROOT}/MTAG"
CONDA_ENV_NAME="mtag_py27"

echo "[INFO] MTAG root directory: ${MTAG_ROOT}"
mkdir -p "${MTAG_ROOT}"

# ---------- Clone MTAG repo ----------
if [[ ! -d "${MTAG_DIR}/mtag" ]]; then
  echo "[INFO] Cloning MTAG repository into ${MTAG_DIR}..."
  mkdir -p "${MTAG_DIR}"
  cd "${MTAG_DIR}"
  git clone https://github.com/omeed-maghzian/mtag.git
else
  echo "[INFO] MTAG repository already present at ${MTAG_DIR}/mtag"
fi

# ---------- Create conda env with Python 2.7 ----------
echo "[INFO] Creating conda environment '${CONDA_ENV_NAME}' (if not existing)..."
if ! conda env list | grep -q "^${CONDA_ENV_NAME}\s"; then
  conda create -y -n "${CONDA_ENV_NAME}" python=2.7
else
  echo "[INFO] Conda environment '${CONDA_ENV_NAME}' already exists."
fi

# ---------- Install Python packages ----------
echo "[INFO] Installing required Python packages into '${CONDA_ENV_NAME}'..."
conda install -y -n "${CONDA_ENV_NAME}" numpy scipy pandas

# Use pip inside the environment for remaining packages
conda run -n "${CONDA_ENV_NAME}" pip install argparse bitarray joblib

# ---------- Test MTAG ----------
echo "[INFO] Testing MTAG help message..."
cd "${MTAG_DIR}"
conda run -n "${CONDA_ENV_NAME}" python mtag/mtag.py -h || {
  echo "[ERROR] MTAG test failed. Check Python / conda setup." >&2
  exit 1
}

echo "[INFO] MTAG installation appears functional."

# ---------- Install CC-GWAS R package ----------
echo "[INFO] Installing CC-GWAS R package from GitHub (wouterpeyrot/CCGWAS)..."

Rscript -e "if (!require('remotes')) install.packages('remotes', repos='https://cloud.r-project.org'); remotes::install_github('wouterpeyrot/CCGWAS')"

echo "[INFO] Verifying CC-GWAS installation..."
Rscript -e "library(CCGWAS); packageVersion('CCGWAS'); sessionInfo()"

echo "[INFO] MTAG + CC-GWAS setup completed successfully."
echo "[INFO] MTAG directory: ${MTAG_DIR}"
echo "[INFO] Conda environment: ${CONDA_ENV_NAME}"

