#!/bin/bash
#

CWD=$(readlink -f $0)
CWD=$(dirname ${CWD})
cd ${CWD}

CONDA_ENV="BIO38"
. ${HOME}/.bashrc
export PATH=${HOME}/local/anaconda3/bin:${HOME}/local/bin:${PATH}

. activate ${CONDA_ENV} || exit 1

PYTHON_MODULE_PATH="/public/dbu/CASP/alphafold/"
export PYTHONPATH="$PYTHONPATH:$PYTHON_MODULE_PATH"



while getopts "i:n:o:" opt; do
  case $opt in
    i) INPUT_FOLDER=$OPTARG ;;
    n) PROTEIN_NAME=$OPTARG ;;
    o) OUTPUT_FOLDER=$OPTARG ;;
    \?) echo "Invalid option: -$OPTARG" >&2
        exit 1 ;;
  esac
done

mkdir -p "$OUTPUT_FOLDER"

# 运行 Python 脚本
python3 "${CWD}/SpliceMSA/Splice-v3.py" -i "$INPUT_FOLDER" -n "$PROTEIN_NAME" -o "$OUTPUT_FOLDER"
