#!/bin/bash

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

cd ${CWD}
echo "command line is:$0 $*"
input_fasta=""
ARGS=$(getopt --options "i:"  -- "$@") || exit
eval "set -- ${ARGS}"
echo "ARGS: $*"
while true; do
  echo "1:$1"
  echo "2:$2"
  case "$1" in
    (-i ) input_fasta="$2"; shift 2;;
    (--) shift 1; break;;
    (*) help 1;;
  esac
done

af2_pred_data_dir="/mnt/data/tanqian/usedeepmsa-dom-casp"
input_fasta="${af2_pred_data_dir}/fasta/${input_fasta}"
echo "input_fasta:${input_fasta}"
af2_params="--use_precomputed_msas --gpu_devices=0,1,3 --multi_gpus --num_predictions_per_model=20"



f=$(basename ${input_fasta})
f=${f%.*}
#for dir in "dMSA" "qMSA" "mMSA/a" "mMSA/b" "mMSA/c"; do
#for dir in "qMSA" "mMSA/a" "mMSA/b" "mMSA/c"; do
for dir in "qMSA"; do
    output_dir="${af2_pred_data_dir}/${f}/${dir}"
    echo "nohup bash -x run_alphafold.sh \
    -m alphafold2:mp \
    -d db \
    -s \
    -o ${output_dir} \
    -p "${af2_params}" \
    ${input_fasta} >> ${af2_pred_data_dir}/${f}.log "

    nohup bash -x run_alphafold.sh \
        -m alphafold2:mp \
        -d db \
        -s \
        -o ${output_dir} \
        -p "${af2_params}" \
        ${input_fasta} \
        >>${af2_pred_data_dir}/${f}.log  

    echo "cleanup:find ${output_dir}/${f} -name "result_model_*.pkl" -exec rm -f {} \;" 
    find ${output_dir}/${f} -name "result_model_*.pkl" -exec rm -f {} \;
    ls -l ${output_dir}/${f}
done



