#!/bin/bash
#

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

cd ${CWD}

help() {
  echo "usage: `basename $0` -d <data_dir> -o <output_dir> -p <msa_params> input_fasta ..."
  exit $1
}

msa_tool_home="/public/dbu"
msa_tool_srv="dbu@idnode17"
msa_tool_sync="rsync -avP -e ssh"
msa_tool_run="ssh ${msa_tool_srv}"

timestamp=$(date +%Y%m%d-%H%M%S-%N)

msa_tool_data_dir="${msa_tool_home}/CASP/data/${timestamp}"

${msa_tool_run} mkdir -p ${msa_tool_data_dir}/fasta

${msa_tool_sync} $* ${msa_tool_srv}:${msa_tool_data_dir}/fasta/


msa_tool_file_list=""
for f in $*; do
  f=$(basename ${f})
  msa_tool_file_list="${msa_tool_file_list} ${msa_tool_data_dir}/fasta/${f}"
done

${msa_tool_run} bash -x ${msa_tool_home}/CASP/run_msa_tool.sh \
    -d ${msa_tool_home}/data/af2_db \
    -o ${msa_tool_data_dir}/msa \
    ${msa_tool_file_list}

echo "msa_tool_run done"

af2_pred_data_dir="${CWD}/data/${timestamp}"

mkdir -p ${af2_pred_data_dir}

${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/msa/* ${af2_pred_data_dir}/

bash -x run_alphafold.sh \
    -m alphafold \
    -d db \
    -o ${af2_pred_data_dir} \
    -p "--use_precomputed_msas" \
    $*
