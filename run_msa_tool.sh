#!/bin/bash
#

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

cd ${CWD}

help() {
  echo "usage: `basename $0` -d <data_dir> -o <output_dir> -p <msa_params> input_fasta ..."
  exit $1
}

export PYTHONPATH=${CWD}/alphafold

########################
# Temp directory
########################
export TMPDIR=tmp
mkdir -p ${TMPDIR}

export JACKHMMER_CPU_NUM=8
export HHBLITS_CPU_NUM=32

CONDA_ENV="af2"
WORKER_NUM=25

docker_image_name="alphafold2"
data_dir=/mnt/data/alphafold/db
output_dir=/mnt/data/alphafold/test
max_template_date="2022-05-01"
msa_params=""
while getopts 'x:n:d:o:m:t:p:h' OPT; do
  case $OPT in
    x) CONDA_ENV="$OPTARG";;
    n) WORKER_NUM="$OPTARG";;
    d) data_dir="$OPTARG";;
    o) output_dir="$OPTARG";;
    p) msa_params="$OPTARG";;
    t) max_template_date="$OPTARG";;
    h) help 0;;
    ?) help 1;;
  esac
done
shift $((OPTIND - 1))

msa_tool_cmd="bash -x $0 -x ${CONDA_ENV} -n ${WORKER_NUM} -d ${data_dir} -o ${output_dir} -t ${max_template_date}"

if [ $# -eq 0 ]; then
  while read input_fasta; do
    ${msa_tool_cmd} -p "${msa_params}" ${input_fasta}
  done
elif [ $# -eq 1 ]; then
  . ${HOME}/.bashrc
  CONDA_HOME=${HOME}/.local/anaconda3
  export PATH=${CONDA_HOME}/bin:${HOME}/.local/bin:${PATH}
  export LD_LIBRARY_PATH=${HOME}/.local/lib:${HOME}/local/lib:${LD_LIBRARY_PATH}
  . activate ${CONDA_ENV} || exit 1

  mkdir -p ${output_dir}

  pid=$(basename $1 | awk -F '.' '{printf("%s", $1); for(i=2;i<NF;++i){printf(".%s", $i);}}');
  find ${output_dir}/${pid}/msas -name "*.sto.gz" -exec gunzip {} \;
  ${CONDA_HOME}/envs/${CONDA_ENV}/bin/python run_msa_tool.py \
        --output_dir=${output_dir} \
        --uniref90_database_path=${data_dir}/uniref90/uniref90.fasta \
        --mgnify_database_path=${data_dir}/mgnify/mgy_clusters.fa \
        --bfd_database_path=${data_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
        --uniref30_database_path=${data_dir}/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
        --template_mmcif_dir=${data_dir}/pdb_mmcif/mmcif_files \
        --obsolete_pdbs_path=${data_dir}/pdb_mmcif/obsolete.dat \
        --pdb70_database_path=${data_dir}/pdb70/pdb70 \
        --pdb_seqres_database_path=${data_dir}/pdb_seqres.txt \
        --uniprot_database_path=${data_dir}/uniprot/uniprot.fasta \
        --max_template_date=${max_template_date} \
        ${msa_params} \
        --fasta_paths $1
else
  if [ $# -lt ${WORKER_NUM} ]; then
    WORKER_NUM=$#
  fi

  echo $@ | \
    awk '{for (i=1; i<=NF; ++i){print $i}}' | \
    ./bin/mapred -m "${msa_tool_cmd} -p \"${msa_params}\"" -c ${WORKER_NUM}
fi
