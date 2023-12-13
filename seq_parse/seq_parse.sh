#!/bin/bash
# 
# This script is used to parse -i <fasta file> to domain fasta files in -o <output_folder>
# use FUpred by zhangyang's lab to predict domain boundary

CWD=$(readlink -f $0)
CWD=$(dirname ${CWD})
cd ${CWD}

CONDA_ENV="BIO38"
. ${HOME}/.bashrc
export PATH=${HOME}/local/anaconda3/bin:${HOME}/local/bin:${PATH}
. activate ${CONDA_ENV} || exit 1

python --version

while getopts "i:o:" opt; do
  case $opt in
    i) INPUT_FASTA=$OPTARG ;;
    o) OUTPUT_FOLDER=$OPTARG ;;
    \?) echo "Invalid option: -$OPTARG" >&2
        exit 1 ;;
  esac
done

mkdir -p "$OUTPUT_FOLDER"

# 运行 seq_parse.py
python3 "${CWD}/seq_parse.py" -i "$INPUT_FASTA" -o "$OUTPUT_FOLDER"