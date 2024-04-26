#!/bin/bash
#

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

cd ${CWD}

timestamp=$(date +%Y%m%d-%H%M%S-%N)
docker_image_name="alphafold2"
output_dir=${CWD}/data
num_predictions_per_model=50

bash ${CWD}/run_casp.sh -m ${docker_image_name} \
	-o ${output_dir} \
	-t "casp-multimer-${timestamp}" \
	--gpu_devices=all \
	--multi_gpus \
	-p "--model_preset=multimer --num_multimer_predictions_per_model=${num_predictions_per_model} --use_precomputed_msas" \
	$*
