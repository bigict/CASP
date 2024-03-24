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
          "  -g, --gpu_devices <gpu_devices>\n" \
          "  -s, --enable_domain_split\n" \
          " --num_predictions_per_model <num_predictions_per_model>\n" \
          "\n" \
          "platforms:\n" \
          "  --msa_tool_home <msa_tool_home>\n" \
          "  --msa_tool_srv <msa_tool_srv>\n" \
          "\n" \
          "workflow:\n" \
          "  [MSA preprocess]* -> MSA search -> [MSA postprocess]* -> Predict"
  exit $1
}

msa_tool_home="/public/dbu"
msa_tool_srv="dbu@idnode17"
msa_tool_sync="rsync -avP -e ssh"



timestamp=$(date +%Y%m%d-%H%M%S-%N)
output_dir=${CWD}/data
model_params=""
gpu_devices="all"
WORKER_NUM=25
enable_domain_split=0
num_predictions_per_model=1 #default

ARGS=$(getopt --options "t:o:p:g:sh" --longoptions "timestamp:,output_dir:,model_params:,gpu_devices:,enable_domain_split,help,msa_tool_home:,msa_tool_srv:,num_predictions_per_model:" -- "$@") || exit
eval "set -- ${ARGS}"
echo "ARGS: $*"
while true; do
  echo "1:$1"
  echo "2:$2"
  case "$1" in
    (-t | --timestamp) timestamp="$2"; shift 2;;
    (-o | --output_dir) output_dir="$2"; shift 2;;
    (-p | --model_params) model_params="$2"; shift 2;;
    (-g | --gpu_devices) gpu_devices="$2"; shift 2;;
    (-h | --help) help 0 ;;
    (-s | --enable_domain_split) enable_domain_split=1; shift 1;;
    (--msa_tool_home) msa_tool_home="$2"; shift 2;;
    (--msa_tool_srv) msa_tool_srv="$2"; shift 2;;
    (--num_predictions_per_model) num_predictions_per_model="$2"; shift 2;;
    (--) shift 1; break;;
    (*) help 1;
  esac
done
echo "timestamp: ${timestamp}"
echo "enable_domain_split: ${enable_domain_split}"
echo "model_params: ${model_params}"
echo "gpu_devices: ${gpu_devices}"
echo "msa_tool_home: ${msa_tool_home}"
echo "msa_tool_srv: ${msa_tool_srv}"
echo "num_predictions_per_model: ${num_predictions_per_model}"
msa_tool_run="ssh ${msa_tool_srv}"


msa_tool_data_dir="${msa_tool_home}/CASP/data/${timestamp}"
${msa_tool_run} mkdir -p ${msa_tool_data_dir}/fasta

#########################
# MSA preprocess modules
#########################

# fasta上传远程
${msa_tool_sync} $* ${msa_tool_srv}:${msa_tool_data_dir}/fasta/

af2_pred_data_dir="${output_dir}/${timestamp}"
mkdir -p ${af2_pred_data_dir}


# 创建msa_tool_file_list文件，里面存放的是各个fasta文件,空格分隔
msa_tool_file_list=""
if [ ${enable_domain_split} -eq 0 ]; then
  for f in $*; do
    f=$(basename ${f})
    msa_tool_file_list="${msa_tool_file_list} ${msa_tool_data_dir}/fasta/${f}"
  done


fi

if [ $# -lt ${WORKER_NUM} ]; then
  WORKER_NUM=$#
fi

if [ ${enable_domain_split} -eq 1 ]; then
  mapred_file_list=""
  for f in $@;do
    f=$(basename ${f})
    mapred_file_list="${mapred_file_list} ${msa_tool_data_dir}/fasta/${f}"
  done
  echo "mapred_file_list: ${mapred_file_list}"
  seq_parse_cmd=" \
  ${msa_tool_run} bash ${msa_tool_home}/CASP/seq_parse/seq_parse.sh \
  -o ${msa_tool_data_dir}/fasta \
  -i"

  echo "${mapred_file_list}" |
  awk '{for (i=1; i<=NF; ++i){print $i}}' | \
  xargs -t -n 1 -P ${WORKER_NUM} -I {} ${seq_parse_cmd} {}
  #./bin/mapred -m "${seq_parse_cmd}" -c ${WORKER_NUM}
  # 远程对input1_fasta input2_fasta ...run seq_parser.py
  # for f in $@; do
  #   f=$(basename ${f})
  #   ${msa_tool_run} bash -x ${msa_tool_home}/CASP/seq_parse/seq_parse.sh -i ${msa_tool_data_dir}/fasta/${f} -o ${msa_tool_data_dir}/fasta
  # done

  echo "parse_domain done"

  # 将远程的${msa_tool_data_dir}/fasta分好的同步到本地
  #parse_data_dir="${af2_pred_data_dir}/parse"
  #mkdir -p ${parse_data_dir}
  ${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/fasta ${af2_pred_data_dir}/

  for f in "${af2_pred_data_dir}"/fasta/* ; do
    f=$(basename ${f})
    #echo "f:${f}"
    msa_tool_file_list="${msa_tool_file_list} ${msa_tool_data_dir}/fasta/${f}"
  done

  echo "msa_tool_file_list: ${msa_tool_file_list}"
fi


echo "start run msa_tool"

msa_params=""
if [ x"${model_params}" == x"" ]; then
  msa_params="-p \"${model_params}\""
fi

${msa_tool_run} bash  ${msa_tool_home}/CASP/run_msa_tool.sh \
    -d ${msa_tool_home}/data/af2_db \
    -o ${msa_tool_data_dir}/msa \
    "${msa_params}" \
    ${msa_tool_file_list}

echo "msa_tool_run done"

#########################
# MSA postprocess modules
#########################
if [ ${enable_domain_split} -eq 1 ]; then
  # run msa_assemble
  for f in $@; do
    f=$(basename ${f})
    f=${f%.*}
    ${msa_tool_run} bash  ${msa_tool_home}/CASP/msa_assemble.sh \
      -i "${msa_tool_data_dir}/msa" \
      -n "${f}" \
      -o "${msa_tool_data_dir}/msa/${f}_assembled/msas"
  done
  echo "msa_assemble done"
fi

${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/fasta ${af2_pred_data_dir}
${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/msa ${af2_pred_data_dir}



#########################
# running alphafold
#########################
if [ x"${model_params}" != x"" ]; then
  model_params="--use_precomputed_msas"
  echo "model_params:${model_params}"
else
  has_use_precomputed_msas=$(echo ${model_params}|grep use_precomputed_msas|wc -l)
  if [ ${has_use_precomputed_msas} -eq 0 ]; then
    model_params="${model_params} --use_precomputed_msas"
  fi
fi
model_params="${model_params} --gpu_devices=${gpu_devices} --num_predictions_per_model=${num_predictions_per_model}"

af2_pred_file_list=""
for f in $*; do
  f=$(basename ${f})
  af2_pred_file_list="${af2_pred_file_list} ${af2_pred_data_dir}/fasta/${f}"
done

if [ ${enable_domain_split} -eq 1 ]; then
  echo "prepare for running alphafold"

  for f in ${af2_pred_file_list}; do 
    f=$(basename ${f})
    echo "f:${f}"
    f=${f%.*}
    #-------------deal with fasta-----------------
    # #把msa_tool_data_dir/fasta中的含"${f}Domain"的fasta文件删除
    # ${msa_tool_run} rm -f ${msa_tool_data_dir}/fasta/${f}Domain*.fasta
    # #copy msa_tool_data_dir/fasta中"${f}.fasta"到msa_tool_data_dir/fasta中，命名为"${f}_assembled.fasta"
    # ${msa_tool_run} cp ${msa_tool_data_dir}/fasta/${f}.fasta ${msa_tool_data_dir}/fasta/${f}_assembled.fasta
    # #同步msa_tool_data_dir/fasta到本地
    # ${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/fasta ${af2_pred_data_dir}
    # #讲${F}_assembled.fasta加入af2_pred_file_list                                                                                                                                                                                                                                                                 
    # af2_pred_file_list="${af2_pred_file_list} ${af2_pred_data_dir}/fasta/${f}_assembled.fasta"


    #把本地的${af2_pred_data_dir}/fasta中"${f}Domain"的fasta文件删除
    rm -f ${af2_pred_data_dir}/fasta/${f}Domain*.fasta
    cp ${af2_pred_data_dir}/fasta/${f}.fasta ${af2_pred_data_dir}/fasta/${f}_assembled.fasta
    #将{f}_assembled.fasta加入af2_pred_file_list                                                                                                                                                                                                                                                                 
    af2_pred_file_list="${af2_pred_file_list} ${af2_pred_data_dir}/fasta/${f}_assembled.fasta"

    echo "af2_pred_file_list:${af2_pred_file_list}"
    #-------------deal with msa-----------------
    rm -rf ${af2_pred_data_dir}/msa/${f}Domain*

    #把${af2_pred_data_dir}/msa/${f}_assembled/msas文件夹中的${f}-assembled.a3m替换同文件目录下的bfd_uniref_hits.a3m
    assembled_dir="${af2_pred_data_dir}/msa/${f}_assembled/msas"
    echo "assembled_dir:${assembled_dir}"
    ls -l "${assembled_dir}"
    mv "${assembled_dir}/${f}-assembled.a3m" "${assembled_dir}/bfd_uniref_hits.a3m"
    ls -l "${assembled_dir}"
    
  done

fi

#delete orginal protein prediction, to speed up
if [ ${enable_domain_split} -eq 1 ]; then
  ls -l ${af2_pred_data_dir}/fasta
  for f in $*; do
    f=$(basename ${f})
    rm -f ${af2_pred_data_dir}/fasta/${f}
  done
  ls -l ${af2_pred_data_dir}/fasta
  #更新af2_pred_list
  af2_pred_file_list=""
  for f in ${af2_pred_data_dir}/fasta/*; do
    f=$(basename ${f})
    af2_pred_file_list="${af2_pred_file_list} ${af2_pred_data_dir}/fasta/${f}"
  done
  echo "to speed up,now af2_pred_file_list:${af2_pred_file_list}"
fi
echo "start run alphafold"
bash -x run_alphafold.sh \
    -m alphafold \
    -d db \
    -o "${af2_pred_data_dir}/msa" \
    -s \
    -p "${model_params}" \
    ${af2_pred_file_list}

# s: single: 让sto生成一个空的；一般是确认我之前搜的msa是否正确