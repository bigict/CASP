#!/bin/bash
#

set -e

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

cd ${CWD}

help() {
  echo -e "usage: `basename $0` [options] input_fasta1 input_fasta2 ...\n" \
          "options:\n" \
          "  -t, --timestamp <timestamp>\n" \
          "  -f, --force\n" \
    "\n" \
    "platforms:\n" \
    "  --msa_tool_home <msa_tool_home>\n" \
    "  --msa_tool_srv <msa_tool_srv>\n"
  exit $1
}

msa_tool_home="/public/dbu"
msa_tool_srv="dbu@inode21"
num_seeds=100
norun_data_pipeline=0
norun_inference=0

timestamp="af3-$(date +%Y%m%d-%H%M%S-%N)"
output_dir=${CWD}/data

ARGS=$(getopt --options "t:fh" --longoptions "timestamp:,msa_tool_home:,msa_tool_srv:,num_seeds:,norun_data_pipeline,norun_inference,help" -- "$@") || exit
eval "set -- ${ARGS}"
echo "ARGS: $*"
while true; do
  case "$1" in
    (-t | --timestamp) timestamp="$2"; shift 2;; 
    (--msa_tool_home) msa_tool_home="$2"; shift 2;; 
    (--msa_tool_srv) msa_tool_srv="$2"; shift 2;; 
    (--num_seeds) num_seeds="$2"; shift 2;; 
    (--norun_data_pipeline) norun_data_pipeline=1; shift 1;; 
    (--norun_inference) norun_inference=1; shift 1;; 
    (-h | --help) help 0 ;;
    (--) shift 1; break;;
    (*) help 1;
  esac
done

msa_tool_sync="rsync -avrzP -e ssh"
msa_tool_run="ssh ${msa_tool_srv}"

msa_tool_data_dir="${msa_tool_home}/alphafold3/data/${timestamp}"
${msa_tool_run} mkdir -p ${msa_tool_data_dir}/fasta

#########################
# MSA preprocess modules
#########################

# fasta上传远程
${msa_tool_sync} $* ${msa_tool_srv}:${msa_tool_data_dir}/fasta/


# 创建msa_tool_file_list文件，里面存放的是各个fasta文件,空格分隔
msa_tool_file_list=""
for f in $*; do
  f=$(basename ${f})
  msa_tool_file_list="${msa_tool_file_list} ${msa_tool_data_dir}/fasta/${f}"
done

if [ ${norun_data_pipeline} -eq 0 ]; then
  echo "start run msa_tool"
  
  ${msa_tool_run} bash -x ${msa_tool_home}/alphafold3/run_pipeline.sh \
      -o ${msa_tool_data_dir}/msa \
      ${msa_tool_file_list}
  
  echo "msa_tool_run done"
fi

af3_pred_data_dir="${output_dir}/${timestamp}"
mkdir -p ${af3_pred_data_dir}

${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/fasta ${af3_pred_data_dir}
${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/msa ${af3_pred_data_dir}

af3_pred_file_list=""
for f in $*; do
  f=$(basename ${f})
  af3_pred_file_list="${af3_pred_file_list} ${af3_pred_data_dir}/fasta/${f}"
done

echo ${af3_pred_file_list}

if [ ${num_seeds} -gt 0 ]; then
    num_seeds_args="--num_seeds=${num_seeds}"
fi

if [ ${norun_inference} -eq 0 ]; then
  # To work around a known XLA issue causing the compilation time to greatly
  # increase, the following environment variable setting XLA flags must be enabled
  # when running AlphaFold 3. Note that if using CUDA capability 7 GPUs, it is
  # necessary to set the following XLA_FLAGS value instead:
  export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
  # (no need to disable gemm in that case as it is not supported for such GPU).
  export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
  # Memory settings used for folding up to 5,120 tokens on A100 80 GB.
  # export XLA_PYTHON_CLIENT_PREALLOCATE=true
  # export XLA_CLIENT_MEM_FRACTION=0.95
  
  PATH=alphafold3/.venv/bin:${PATH} python alphafold3/run_alphafold.py \
    --norun_data_pipeline \
    --jax_compilation_cache_dir ${CWD}/.cache/jax_compilation \
    --model_dir ${HOME}/public_databases/models \
    --input_dir ${af3_pred_data_dir}/msa \
    --output_dir ${af3_pred_data_dir}/msa \
    ${num_seeds_args}
fi
