#!/bin/bash
#

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

cd ${CWD}

help() {
  echo "usage: `basename $0` -d <data_dir> -o <output_dir> -p <msa_params> input_fasta ..."
  exit $1
}

docker_image_name="alphafold2"
data_dir=/mnt/data/alphafold/db
output_dir=/mnt/data/alphafold/test
max_template_date="2021-05-14"
msa_params=""
while getopts 'd:o:m:t:p:h' OPT; do
  case $OPT in
    d) data_dir="$OPTARG";;
    o) output_dir="$OPTARG";;
    p) msa_params="$OPTARG";;
    t) max_template_date="$OPTARG";;
    h) help 0;;
    ?) help 1;;
  esac
done
shift $((OPTIND - 1))

mkdir -p ${output_dir}
for f in $*; do
  echo "Start predict (${f}) `date`"
  python run_msa_tool.py \
      --output_dir=${output_dir} \
      --uniref90_database_path=${data_dir}/uniref90/uniref90.fasta \
      --mgnify_database_path=${data_dir}/mgnify/mgy_clusters.fa \
      --bfd_database_path=${data_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
      --uniclust30_database_path=${data_dir}/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
      --template_mmcif_dir=${data_dir}/pdb_mmcif/mmcif_files \
      --obsolete_pdbs_path=${data_dir}/pdb_mmcif/obsolete.dat \
      --pdb70_database_path=${data_dir}/pdb70/pdb70 \
      --max_template_date=${max_template_date} \
      --fasta_paths=${f} \
      ${msa_params}
  echo "Finished (${f}) `date`"
done
