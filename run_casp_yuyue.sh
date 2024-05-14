if [ ${force_predict} -ne 0 -o ! -e ${output_dir}/${timestamp}/DONE ]; then
  bash -x ${CWD}/run_casp.sh -m ${docker_image_name} \
        -o ${output_dir} \
        -t "${timestamp}" \
        --gpu_devices=all \
        --multi_gpus \
        -p "--model_preset=monomer_ptm --num_predictions_per_model=${num_predictions_per_model} --use_precomputed_msas" \
        $*
fi

for fasta_file in $*; do
   pid=$(basename ${fasta_file} | awk -F '.fa|.fasta' '{print $1}')
   python massivefold/massivefold/massivefold_plots.py \
        --chosen_plots coverage,CF_plddt,DM_plddt_PAE,CF_PAEs,CF_plddts,score_distribution \
        --input_path ${output_dir}/${timestamp}/msa/${pid}
done
