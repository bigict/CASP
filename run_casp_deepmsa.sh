#!/bin/bash
#

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

cd ${CWD}
echo "the whole command line is:$0 $*"

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
msa_tool_srv="dbu@idnode19"
msa_tool_sync="rsync -avP -e ssh"



timestamp=$(date +%Y%m%d-%H%M%S-%N)
#output_dir=${CWD}/data
output_dir="/mnt/data/tanqian"
model_params=""
gpu_devices="all"
WORKER_NUM=25
enable_domain_split=0
num_predictions_per_model=1 #default
clean_model_result_pkl=0
multi_gpus=0
use_deepmsa2=0
skip_af2_msa=0
ARGS=$(getopt --options "t:o:p:g:shc" --longoptions "timestamp:,output_dir:,model_params:,gpu_devices:,enable_domain_split,help,msa_tool_home:,msa_tool_srv:,num_predictions_per_model:,multi_gpus,skip_af2_msa,use_deepmsa2" -- "$@") || exit
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
    (-c ) clean_model_result_pkl=1; shift 1;;
    (-s | --enable_domain_split) enable_domain_split=1; shift 1;;
    (--msa_tool_home) msa_tool_home="$2"; shift 2;;
    (--msa_tool_srv) msa_tool_srv="$2"; shift 2;;
    (--num_predictions_per_model) num_predictions_per_model="$2"; shift 2;;
    (--multi_gpus) multi_gpus=1; shift 1;;
    (--skip_af2_msa) skip_af2_msa=1;shift 1;;
    (--use_deepmsa2) use_deepmsa2=1; shift 1;;
    (--num_predictions_per_model) num_predictions_per_model="$2"; shift 2;;
    (--) shift 1; break;;
    (*) help 1;
  esac
done

echo "timestamp: ${timestamp}"
echo "output_dir:${output_dir}"
echo "enable_domain_split: ${enable_domain_split}"
echo "model_params: ${model_params}"
echo "gpu_devices: ${gpu_devices}"
echo "msa_tool_home: ${msa_tool_home}"
echo "msa_tool_srv: ${msa_tool_srv}"
echo "num_predictions_per_model: ${num_predictions_per_model}"
echo "clean_model_result_pkl: ${clean_model_result_pkl}"
echo "multi_gpus: ${multi_gpus}"
echo "skip_af2_msa: ${skip_af2_msa}"
echo "use_deepmsa2: ${use_deepmsa2}"
msa_tool_run="ssh ${msa_tool_srv}"


msa_tool_data_dir="${msa_tool_home}/CASP_multimer/data/${timestamp}"
${msa_tool_run} mkdir -p ${msa_tool_data_dir}/fasta
af2_pred_data_dir="${output_dir}/${timestamp}"
mkdir -p ${af2_pred_data_dir}

#########################
# MSA preprocess modules
#########################

# fasta上传远程
${msa_tool_sync} $* ${msa_tool_srv}:${msa_tool_data_dir}/fasta/



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
  ${msa_tool_run} bash ${msa_tool_home}/CASP_multimer/seq_parse/seq_parse.sh \
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

if [ ${skip_af2_msa} -eq 0 ]; then
  msa_params=""
  if [ x"${model_params}" != x"" ]; then
    msa_params="-p \"${model_params}\""
  fi

  echo "start run msa_tool"
  ${msa_tool_run} bash  ${msa_tool_home}/CASP_multimer/run_msa_tool.sh \
    -d ${msa_tool_home}/data/af2_db \
    -o ${msa_tool_data_dir}/msa \
    "${msa_params}" \
    ${msa_tool_file_list}
  echo "run_msa_tool_run done"
else
  echo "skip af2 searching MSA"
fi

${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/fasta ${af2_pred_data_dir}
${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/msa ${af2_pred_data_dir}


#########################
# prepare for af2:params
#########################
if [ x"${model_params}" == x"" ]; then
  model_params="--use_precomputed_msas"
  echo "model_params:${model_params}"
else
  has_use_precomputed_msas=$(echo ${model_params}|grep use_precomputed_msas|wc -l)
  if [ ${has_use_precomputed_msas} -eq 0 ]; then
    model_params="${model_params} --use_precomputed_msas"
  fi
fi

if [ ${multi_gpus} -eq 1 ]; then
  model_params="${model_params} --multi_gpus"
fi
model_params="${model_params} --gpu_devices=${gpu_devices} --num_predictions_per_model=${num_predictions_per_model}"
echo "model_params:${model_params}"



#########################
# MSA postprocess modules
#########################

if [ ${use_deepmsa2} -eq 1 ]; then  
  echo "start run deepmsa2"
  #统计${msa_tool_file_list}中含有字符串“Domain”的fasta文件数
  deepmsa2_work_num=$(echo ${msa_tool_file_list}|grep -o "Domain"|wc -l)
  deepmsa2_file_list=""
  for f in ${msa_tool_file_list}; do
    if [ $(echo ${f}|grep -o "Domain"|wc -l) -eq 1 ]; then
      deepmsa2_file_list="${deepmsa2_file_list} ${f}"
    fi
  done
  echo "deepmsa2_file_list:${deepmsa2_file_list}"
  echo "${msa_tool_run} bash  ${msa_tool_home}/CASP_multimer/run_deep_msa2.sh \
    -o ${msa_tool_data_dir}/msa \
    -n ${deepmsa2_work_num} \
    ${deepmsa2_file_list}"

  # ${msa_tool_run} bash ${msa_tool_home}/CASP_multimer/run_deep_msa2.sh \
  #   -o ${msa_tool_data_dir}/msa \
  #   -n ${deepmsa2_work_num} \
  #   ${deepmsa2_file_list}
  echo "deepmsa2 search msa done"
  
  #convert deepmsa2 format to af2 format,running af2
  for f in ${deepmsa2_file_list}; do
    fasta_name=$(basename ${f})
    fasta_name=${fasta_name%.*}
    echo "${msa_tool_run} bash  ${msa_tool_home}/CASP_multimer/deepmsa2af.sh \
      -o "${msa_tool_data_dir}/msa/${fasta_name}" \
      --hhr_file "/public/dbu/CASP_multimer/data/424-native/msa/${fasta_name}/msas/pdb_hits.hhr" \
      ${f}"

    ${msa_tool_run} bash  ${msa_tool_home}/CASP/deepmsa2af.sh \
      -o "${msa_tool_data_dir}/msa/${fasta_name}" \
      --hhr_file "/public/dbu/CASP_multimer/data/424-native/msa/${fasta_name}/msas/pdb_hits.hhr" \
      ${f}
    echo "${f} deepmsa2af done"
  done

  ${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/fasta ${af2_pred_data_dir}
  ${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/msa ${af2_pred_data_dir}
  

  for f in ${deepmsa2_file_list}; do
    f=$(basename ${f})
    fasta_file=${af2_pred_data_dir}/fasta/${f}
    f=${f%.*} #f="T1170Domain1"
    declare -A tmscore_dict #关联数组
    for dir in "dMSA" "qMSA" "mMSA/a" "mMSA/b" "mMSA/c"; do
      output_dir="${af2_pred_data_dir}/msa/${f}/${dir}"
      echo "bash -x run_alphafold.sh \
          -m alphafold2:latest \
          -d db\
          -s \
          -o ${output_dir} \
          -p "${model_params}" \
          ${fasta_file} \
          >${output_dir}/${f}.log"

      bash -x run_alphafold.sh \
          -m alphafold2:latest \
          -s \
          -d db\
          -o ${output_dir} \
          -p "${model_params}" \
          ${fasta_file} \
          >${output_dir}/${f}.log

      echo "cleanup:find ${output_dir}/${f} -name "result_model_*.pkl" -exec rm -f {} \;" 
      find ${output_dir}/${f} -name "result_model_*.pkl" -exec rm -f {} \;

      pushd ${output_dir}/${f}
      log_file="${output_dir}/find_max_tmscore.log"
      > ${log_file}
      for pdb in $(ls ranked_*.pdb); do
          pdb=$(readlink -f ${pdb})
          python /mnt/data/tanqian/tmscore.py   --gdt -e /home/tanqian/DeepAlign/DeepScore  -v ${pdb} /mnt/data/tanqian/data/"${f}".pdb >>  ${log_file}
      done
      
      # 使用 grep 命令找出所有 tmscore 值，并排序找出最大值
      max_tmscore=$(grep -o 'tmscore=[0-9.]*' "$log_file" | cut -d= -f2 | sort -nr | head -n1)
      echo "${dir}'s max_tmscore=$max_tmscore"
      popd
      #find max tmscore
      tmscore_dict[$dir]=${max_tmscore}
    done

    #find max tmscore in tmscore_dict
    max_tmscore=0
    max_dir=""
    for dir in ${!tmscore_dict[@]}; do
      tmscore=${tmscore_dict[$dir]}
      if [ $(echo "${tmscore} > ${max_tmscore}"|bc) -eq 1 ]; then
        max_tmscore=${tmscore}
        max_dir=${dir}
      fi
    done
    echo "${f}'s max_tmscore:${max_tmscore} in ${max_dir}"
    mkdir -p ${af2_pred_data_dir}/deepmsa/${f}/msas
    cp -r ${af2_pred_data_dir}/msa/${f}/${max_dir}/${f}/msas ${af2_pred_data_dir}/deepmsa/${f}/  #现在${af2_pred_data_dir}/deepmsa/${f}有msas文件夹

  done
  #delete ${af2_pred_data_dir}/msa,change /deepmsa to /msa一些替身文学hhh
  mv ${af2_pred_data_dir}/msa ${af2_pred_data_dir}/trash-msa
  mv ${af2_pred_data_dir}/deepmsa ${af2_pred_data_dir}/msa
  #将msa_tool_data_dir/msa改名为trash-msa
  ${msa_tool_run} mv ${msa_tool_data_dir}/msa ${msa_tool_data_dir}/trash-msa
  ${msa_tool_sync} ${af2_pred_data_dir}/msa ${msa_tool_srv}:${msa_tool_data_dir}/
fi



if [ ${enable_domain_split} -eq 1 ]; then
  # run msa_assemble
  for f in $@; do
    f=$(basename ${f})
    f=${f%.*}
    echo "${msa_tool_run} bash  ${msa_tool_home}/CASP_multimer/msa_assemble.sh \
      -i "${msa_tool_data_dir}/msa" \
      -n "${f}" \
      -o "${msa_tool_data_dir}/msa/${f}_assembled/msas" \
      --input_fasta "${msa_tool_data_dir}/fasta/${f}.fasta"    "  

    ${msa_tool_run} bash  ${msa_tool_home}/CASP_multimer/msa_assemble.sh \
      -i "${msa_tool_data_dir}/msa" \
      -n "${f}" \
      -o "${msa_tool_data_dir}/msa/${f}_assembled/msas" \
      --input_fasta "${msa_tool_data_dir}/fasta/${f}.fasta"

    
    ${msa_tool_run} cp /public/dbu/CASP_multimer/data/424-native/msa/${f}/msas/pdb_hits.hhr ${msa_tool_data_dir}/msa/${f}_assembled/msas
  done
  echo "msa_assemble done"
fi

${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/fasta ${af2_pred_data_dir}
${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/msa ${af2_pred_data_dir}



#########################
# running alphafold
#########################
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
    #cp ${af2_pred_data_dir}/fasta/${f}.fasta ${af2_pred_data_dir}/fasta/${f}_assembled.fasta
    #将{f}_assembled.fasta加入af2_pred_file_list                                                                                                                                                                                                                                                                 
    #af2_pred_file_list="${af2_pred_file_list} ${af2_pred_data_dir}/fasta/${f}_assembled.fasta"

    echo "af2_pred_file_list:${af2_pred_file_list}"
    #-------------deal with msa-----------------
    rm -rf ${af2_pred_data_dir}/msa/${f}Domain*

    #把${af2_pred_data_dir}/msa/${f}_assembled/msas文件夹中的${f}-assembled.a3m替换同文件目录下的bfd_uniref_hits.a3m
    assembled_dir="${af2_pred_data_dir}/msa/${f}_assembled/msas"
    echo "assembled_dir:${assembled_dir}"
    #ls -l "${assembled_dir}"
    mv "${assembled_dir}/${f}-assembled.a3m" "${assembled_dir}/bfd_uniref_hits.a3m"
    #ls -l "${assembled_dir}"
    mkdir -p ${af2_pred_data_dir}/msa/${f}/msas
    mv ${assembled_dir} ${af2_pred_data_dir}/msa/${f}
    echo "mv ${assembled_dir} ${af2_pred_data_dir}/msa/${f}"
    ls -l ${af2_pred_data_dir}/msa/${f}/msas
  done

fi


echo "start run alphafold"
#db=/mnt/data/alphafold   #node10
echo "bash -x run_alphafold.sh \
    -m alphafold2:latest \
    -o "${af2_pred_data_dir}/msa" \
    -s \
    -d db\
    -p "${model_params}" \
    ${af2_pred_file_list}"

bash -x run_alphafold.sh \
    -m alphafold2:latest \
    -o "${af2_pred_data_dir}/msa" \
    -s \
    -d db\
    -p "${model_params}" \
    ${af2_pred_file_list}
# #db=/mnt/data/alphafold   #node10

if [ ${enable_domain_split} -eq 1 ]; then
  for f  in ${af2_pred_file_list}; do
    f=$(basename ${f})
    f=${f%.*}
    pushd ${af2_pred_data_dir}/msa/${f}
    log_file="${af2_pred_data_dir}/${f}_ass_find_max_tmscore.log"
    > ${log_file}
    for pdb in $(ls ranked_*.pdb); do
      pdb=$(readlink -f ${pdb})
      echo "python /mnt/data/tanqian/tmscore.py   --gdt -e /home/tanqian/DeepAlign/DeepScore  -v ${pdb} /mnt/data/tanqian/data/"${f/_assembled/}".pdb >>  ${log_file}"
      python /mnt/data/tanqian/tmscore.py   --gdt -e /home/tanqian/DeepAlign/DeepScore  -v ${pdb} /mnt/data/tanqian/data/"${f/_assembled/}".pdb >>  ${log_file}   
    done
  done
fi

# s: single: 让sto生成一个空的；一般是确认我之前搜的msa是否正确
if [ ${clean_model_result_pkl} -eq 1 ]; then
  for f in ${af2_pred_file_list}; do
    f=$(basename ${f})
    f=${f%.*}
    echo "cleanup:find ${af2_pred_data_dir}/msa/${f} -name "result_model_*.pkl" -exec rm -f {} \;" 
    find ${af2_pred_data_dir}/msa/${f} -name "result_model_*.pkl" -exec rm -f {} \;
  done
fi