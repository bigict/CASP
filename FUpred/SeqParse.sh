#!/bin/bash
#

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

# 运行 SeqParse.py
python3 "${CWD}/SeqParse.py" -i "$INPUT_FASTA" -o "$OUTPUT_FOLDER"