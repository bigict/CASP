#!/bin/bash
#
# This script is used to assemble msa results from msa_tool
# assembling strategy: 
#     1. use alphafold.data.parser's convert_stockholm_to_a3m
#     2. assemble msas vertically from 3 databases: uniref90, uniref50, bfd within domain
#     3. assemble msas horizontally with domains

help(){
  echo "Usage: $0 -i <input_folder> -n <protein_name> -o <output_folder>"
  echo "Options:"
  echo "  -i <input_folder>    input folder of msa results from msa_tool"
  echo "  -n <protein_name>    protein name,eg:T1125"
  echo "  -o <output_folder>   output folder of assembled msa"
  exit $1
}
CWD=$(readlink -f $0)
CWD=$(dirname ${CWD})
cd ${CWD}

CONDA_ENV="BIO38"
. ${HOME}/.bashrc
export PATH=${HOME}/local/anaconda3/bin:${HOME}/local/bin:${PATH}

. activate ${CONDA_ENV} || exit 1

PYTHON_MODULE_PATH="/public/dbu/CASP/alphafold/"
export PYTHONPATH="$PYTHONPATH:$PYTHON_MODULE_PATH"



while getopts "i:n:o:h" opt; do
  case $opt in
    i) INPUT_FOLDER=$OPTARG ;;
    n) PROTEIN_NAME=$OPTARG ;;
    o) OUTPUT_FOLDER=$OPTARG ;;
    h) help 0 ;;
    ?) help 1 ;;
  esac
done
shift $((OPTIND - 1))

mkdir -p "$OUTPUT_FOLDER"

# 运行 Python 脚本
python3 "${CWD}/msa_assemble.py" -i "$INPUT_FOLDER" -n "$PROTEIN_NAME" -o "$OUTPUT_FOLDER"
