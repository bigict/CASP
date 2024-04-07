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
timestamp=$(date +%Y%m%d-%H%M%S-%N)
#output_dir=${CWD}/data
output_dir="/mnt/data/tanqian"
model_params=""
gpu_devices="all"
WORKER_NUM=25
enable_domain_split=0
num_predictions_per_model=1 #default
clean_model_rsult_pkl=0
ARGS=$(getopt --options "t:o:p:g:shc" --longoptions "timestamp:,output_dir:,model_params:,gpu_devices:,enable_domain_split,help,msa_tool_home:,msa_tool_srv:,num_predictions_per_model:" -- "$@") || exit
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
    (-c) clean_model_rsult_pkl=1; shift 1;;
    (-s | --enable_domain_split) enable_domain_split=1; shift 1;;
    (--msa_tool_home) msa_tool_home="$2"; shift 2;;
    (--msa_tool_srv) msa_tool_srv="$2"; shift 2;;
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
echo "clean_model_rsult_pkl: ${clean_model_rsult_pkl}"
msa_tool_run="ssh ${msa_tool_srv}"
af2_pred_data_dir="${output_dir}/${timestamp}"
mkdir -p ${af2_pred_data_dir}

#########################
# running alphafold
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
model_params="${model_params} --gpu_devices=${gpu_devices} --num_predictions_per_model=${num_predictions_per_model}"
echo "model_params:${model_params}"
af2_pred_file_list=""
for f in $*; do
  f=$(basename ${f})
  af2_pred_file_list="${af2_pred_file_list} ${af2_pred_data_dir}/fasta/${f}"
done
echo "start run alphafold"

db=/mnt/data/alphafold
bash -x run_alphafold.sh \
    -m alphafold2:tanqian8 \
    -d ${db} \
    -o "${af2_pred_data_dir}/msa" \
    -s \
    -p "${model_params}" \
    ${af2_pred_file_list}

if [ ${clean_model_rsult_pkl} -eq 1 ]; then
  for f in ${af2_pred_file_list}; do
    f=$(basename ${f})
    f=${f%.*}
    echo "cleanup:find ${af2_pred_data_dir}/msa/${f} -name "result_model_*.pkl" -exec rm -f {} \;" 
    find ${af2_pred_data_dir}/msa/${f} -name "result_model_*.pkl" -exec rm -f {} \;
    ls -l ${af2_pred_data_dir}/msa/${f}
  done
fi