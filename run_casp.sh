#!/bin/bash
#

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

cd ${CWD}

help() {
  echo -e "usage: `basename $0` [options] input_fasta1 input_fasta2 ...\n" \
          "options:\n" \
          "  -t, --timestamp <timestamp>\n" \
          "  -o, --output_dir <output_dir>\n" \
          "  -p, --model_params <model_params>\n" \
          "  -s, --enable_domain_split\n" \
          "\n" \
          "workflow:\n" \
          "  [MSA preprocess]* -> MSA search -> [MSA postprocess]* -> Predict"
  exit $1
}

msa_tool_home="/public/dbu"
msa_tool_srv="dbu@idnode17"
msa_tool_sync="rsync -avP -e ssh"
msa_tool_run="ssh ${msa_tool_srv}"


timestamp=$(date +%Y%m%d-%H%M%S-%N)
output_dir=${CWD}/data
model_params=""

enable_domain_split=0

ARGS=$(getopt --options "t:o:p:sh" --longoptions "timestamp:,output_dir:,model_params:,enable_domain_split,help" -- "$@") || exit
eval "set -- ${ARGS}"

while true; do
  case "$1" in
    (-t | --timestamp) timestamp="$2" shift 2;;
    (-o | --output_dir) output_dir="$2" shift 2;;
    (-p | --model_params) model_params="-p \"$2\"" shift 2;;
    (-h | --help) help 0;;
    (-s | --enable_domain_split) enable_domain_split=1 shift 1;;
    (--) shift 1; break;;
    (*) help 1;
  esac
done

msa_tool_data_dir="${msa_tool_home}/CASP/data/${timestamp}"
${msa_tool_run} mkdir -p ${msa_tool_data_dir}/fasta

#########################
# MSA preprocess modules
#########################

# fasta上传远程
${msa_tool_sync} $* ${msa_tool_srv}:${msa_tool_data_dir}/fasta/

# 创建msa_tool_file_list文件，里面存放的是各个fasta文件,空格分隔
msa_tool_file_list=""
if [ ${enable_domain_split} -eq 0]; then
  for f in $*; do
  f=$(basename ${f})
  msa_tool_file_list="${msa_tool_file_list} ${msa_tool_data_dir}/fasta/${f}"
done
fi


af2_pred_data_dir="${output_dir}/${timestamp}"
mkdir -p ${af2_pred_data_dir}

if [ ${enable_domain_split} -ne 0 ]; then
  # 远程对input1_fasta input2_fasta ...run seq_parser.py
  for f in $@; do
    f=$(basename ${f})
    ${msa_tool_run} bash -x ${msa_tool_home}/CASP/seq_parse/seq_parse.sh -i ${msa_tool_data_dir}/fasta/${f} -o ${msa_tool_data_dir}/fasta
  done

  echo "parse_domain done"

  # 将远程的${msa_tool_data_dir}/fasta分好的同步到本地
  parse_data_dir="${af2_pred_data_dir}/parse"
  mkdir -p ${parse_data_dir}
  ${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/fasta ${parse_data_dir}/

  for f in "${parse_data_dir}"/fasta/* ; do
    f=$(basename ${f})
    #echo "f:${f}"
    msa_tool_file_list="${msa_tool_file_list} ${msa_tool_data_dir}/fasta/${f}"
  done

  echo "msa_tool_file_list: ${msa_tool_file_list}"
fi


echo "start run msa_tool"

${msa_tool_run} bash -x ${msa_tool_home}/CASP/run_msa_tool.sh \
    -d ${msa_tool_home}/data/af2_db \
    -o ${msa_tool_data_dir}/msa \
    ${model_params} \
    ${msa_tool_file_list}

echo "msa_tool_run done"

#########################
# MSA postprocess modules
#########################
if [ ${enable_domain_split} -ne 0 ]; then
  # check whether exists msa ,use msa_tool_file_list
  for f in $@; do
    f=$(basename ${f})
    f=${f%.*}
    ${msa_tool_run} bash -x ${msa_tool_home}/CASP/msa_assemble.sh \
      -i "${msa_tool_data_dir}/msa" \
      -n "${f}" \
      -o "${msa_tool_data_dir}/msa/${f}_assembled"
  done
fi

${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/msa ${af2_pred_data_dir}/

if [ x"${model_params}" == x"" ]; then
  model_params="-p --use_precomputed_msas"
else
  has_use_precomputed_msas=$(echo ${model_params}|grep use_precomputed_msas|wc -l)
  if [ ${has_use_precomputed_msas} -eq 0 ]; then
    model_params="${model_params} --use_precomputed_msas"
  fi
fi

bash -x run_alphafold.sh \
     -m alphafold \
     -d db \
     -o ${af2_pred_data_dir} \
     ${model_params} \
     ${msa_tool_file_list}
