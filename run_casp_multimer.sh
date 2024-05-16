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
	  "  --stoich <stoichiometry> *required*" \
	  "  --renumber_rna_start <renumber_rna_start>" \
	  "  --renumber_dna_start <renumber_dna_start>" \
          "  --force_predict\n" \
          "  --force_plot\n" \
          "  --force_email\n" \
	  "\n" \
	  "platforms:\n" \
	  "  --msa_tool_home <msa_tool_home>\n" \
	  "  --msa_tool_srv <msa_tool_srv>\n"
  exit $1
}

msa_tool_home="/public/dbu"
msa_tool_srv="dbu@idnode18"
# msa_tool_srv="dbu@idnode17"

timestamp="casp-multimer-$(date +%Y%m%d-%H%M%S-%N)"
force_predict=0
force_plot=0
force_email=0
docker_image_name="alphafold2"
output_dir=${CWD}/data
num_predictions_per_model=${num_predictions_per_model:-"50"}
model_to_relax_topk=${model_to_relax_topk:-"5"}
stoichiometry=""
renumber_rna_start=0
renumber_dna_start=-1

ARGS=$(getopt --options "t:h" --longoptions "timestamp:,msa_tool_home:,msa_tool_srv:,stoich:,renumber_rna_start:,renumber_dna_start:,force_predict,force_plot,force_email,help" -- "$@") || exit
eval "set -- ${ARGS}"
echo "ARGS: $*"
while true; do
  echo "1:$1"
  echo "2:$2"
  case "$1" in
    (-t | --timestamp) timestamp="$2"; shift 2;;
    (--stoich) stoichiometry="$2"; shift 2;;
    (--renumber_rna_start) renumber_rna_start="$2"; shift 2;;
    (--renumber_dna_start) renumber_dna_start="$2"; shift 2;;
    (--force_predict) force_predict=1; shift 1;;
    (--force_plot) force_plot=1; shift 1;;
    (--force_email) force_email=1; shift 1;;
    (--msa_tool_home) msa_tool_home="$2"; shift 2;;
    (--msa_tool_srv) msa_tool_srv="$2"; shift 2;;
    (-h | --help) help 0 ;;
    (--) shift 1; break;;
    (*) help 1;
  esac
done

if [ -z "${stoichiometry}" ]; then
  help 1;
fi

if [ ${force_predict} -ne 0 -o ! -e ${output_dir}/${timestamp}/DONE ]; then
  export MODEL_TO_RELAX_TOPK=${model_to_relax_topk}
  bash -x ${CWD}/run_casp.sh -m ${docker_image_name} \
  	-o ${output_dir} \
  	-t "${timestamp}" \
  	--gpu_devices=all \
  	--multi_gpus \
	--msa_tool_srv ${msa_tool_srv} \
  	-p "--model_preset=multimer --num_multimer_predictions_per_model=${num_predictions_per_model} --use_precomputed_msas" \
  	$*
fi

if [ ${force_plot} -ne 0 -o ! -e ${output_dir}/${timestamp}/DONE ]; then
  for fasta_file in $*; do
    pid=$(basename ${fasta_file} | awk -F '.fa|.fasta' '{print $1}')
    python massivefold/massivefold/massivefold_plots.py \
    	--chosen_plots coverage,CF_plddt,DM_plddt_PAE,CF_PAEs,CF_plddts,score_distribution \
    	--input_path ${output_dir}/${timestamp}/msa/${pid}
  done
fi

if [ ${force_email} -ne 0 -o ! -e ${output_dir}/${timestamp}/DONE ]; then
  for fasta_file in $*; do
    pid=$(basename ${fasta_file} | awk -F '.fa|.fasta' '{print $1}')
    mkdir -p ${output_dir}/${timestamp}/msa/${pid}/email
    for ((idx=0; idx<5; ++idx)); do
      model=$(expr ${idx} + 1)
      python submit_casp.py --dry_run -v \
	  --target=${pid} \
	  --model=${model} \
	  --stoich=${stoichiometry} \
	  --renumber_chain \
	  --renumber_rna_start=${renumber_rna_start} \
	  --renumber_dna_start=${renumber_dna_start} \
	  --output_file=${output_dir}/${timestamp}/msa/${pid}/email/model_${model}.pdb \
	  ${output_dir}/${timestamp}/msa/${pid}/ranked_${idx}.cif
    done
  done
fi

touch ${output_dir}/${timestamp}/DONE
