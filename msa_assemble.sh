#!/bin/bash
#this script is used to assemble msa results from msa_tool
#assembling strategy: 
#1. use alphafold.data.parser's convert_stockholm_to_a3m
#2. assemble msas vertically from 3 databases: uniref90, uniref50, bfd within domain
#3. assemble msas horizontally with domains

help(){
  echo "Usage: $0 -i <input_folder> -n <protein_name> -o <output_folder> --input_fasta <input_fasta>"
  echo "Options:"
  echo "  -i <input_folder>    input folder of msa results from msa_tool"
  echo "  -n <protein_name>    protein name,eg:T1125"
  echo "  -o <output_folder>   output folder of assembled msa"
  echo "  --input_fasta <input_fasta>   input fasta file of protein sequence"
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

ARGS=$(getopt --options "i:o:n:h" --longoptions "input_fasta:" -- "$@") || exit
eval "set -- ${ARGS}"
echo "ARGS: $*"
while true; do
  echo "1:$1"
  echo "2:$2"
  case "$1" in
    (-i ) INPUT_FOLDER="$2"; shift 2;;
    (-o ) OUTPUT_FOLDER="$2"; shift 2;;
    (-n ) PROTEIN_NAME="$2"; shift 2;;
    (--input_fasta) INPUT_FASTA="$2"; shift 2;;
    (-h | --help) help 0 ;;
    (--) shift 1; break;;
    (*) help 1;;
  esac
done


mkdir -p "$OUTPUT_FOLDER"

# 运行 Python 脚本
echo python3 "${CWD}/msa_assemble.py" -i "$INPUT_FOLDER" -n "$PROTEIN_NAME" -o "$OUTPUT_FOLDER --input_fasta $INPUT_FASTA"

python3 "${CWD}/msa_assemble.py" -i "$INPUT_FOLDER" -n "$PROTEIN_NAME" -o "$OUTPUT_FOLDER --input_fasta $INPUT_FASTA"
