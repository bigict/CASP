#!/bin/bash
#

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

cd ${CWD}

help() {
  echo "usage: `basename $0` -d <data_dir> -o <output_dir> -p <af2_params> input_fasta ..."
  exit $1
}

docker_image_name="alphafold2"
data_dir=/mnt/data/alphafold/db
output_dir=/mnt/data/alphafold/test
max_template_date="2021-05-14"
af2_params=""
while getopts 'd:o:m:t:p:h' OPT; do
  case $OPT in
    m) docker_image_name="$OPTARG";;
    d) data_dir="$OPTARG";;
    o) output_dir="$OPTARG";;
    p) af2_params="$OPTARG";;
    t) max_template_date="$OPTARG";;
    h) help 0;;
    ?) help 1;;
  esac
done
shift $((OPTIND - 1))

mkdir -p ${output_dir}
for f in $*; do
  echo "Start predict (${f}) `date`"
  python alphafold/docker/run_docker.py \
      --docker_image_name=${docker_image_name} \
      --data_dir=${data_dir} \
      --max_template_date=${max_template_date} \
      --output_dir=${output_dir} \
      --fasta_paths=${f} \
      ${af2_params}
  echo "Finished (${f}) `date`"
done
