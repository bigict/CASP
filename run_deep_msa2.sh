#!/bin/bash
#

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

cd ${CWD}

help() {
  echo "usage: `basename $0` -o <output_dir> -n <worker_num> input_fasta ..."
  exit $1
}

export HHLIB=${CWD}/msa-suite

CONDA_ENV="af2"
WORKER_NUM=5

output_dir=/mnt/data/alphafold/test
while getopts 'x:n:o:h' OPT; do
  case $OPT in
    x) CONDA_ENV="$OPTARG";;
    n) WORKER_NUM="$OPTARG";;
    o) output_dir="$OPTARG";;
    h) help 0;;
    ?) help 1;;
  esac
done
shift $((OPTIND - 1))


. ${HOME}/.bashrc
CONDA_HOME=${HOME}/.local/anaconda3
export PATH=${CONDA_HOME}/bin:${HOME}/.local/bin:${PATH}
export LD_LIBRARY_PATH=${HOME}/.local/lib:${HOME}/local/lib:${LD_LIBRARY_PATH}
. activate ${CONDA_ENV} || exit 1

mkdir -p ${output_dir}

deep_msa2_cmd="bash $0 -x ${CONDA_ENV} -o ${output_dir}"
if [ $# -eq 0 ]; then
  while read input_fasta; do
    ${deep_msa2_cmd} -n 1 ${input_fasta}
  done
elif [ $# -eq 1 ]; then
  db_dir=${CWD}/msa-suite/db
  export qmsa_bfddb="${db_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt,${db_dir}/uniref30/UniRef30_2022_02"
  pid=$(basename $1 | awk -F '.' '{printf("%s", $1); for(i=2;i<NF;++i){printf(".%s", $i);}}');
  bash -x msa-suite/examples/deepMSA2.sh -o ${output_dir}/${pid} -i $1 -vk
  for msa in "dMSA" "qMSA" "mMSA/a" "mMSA/b" "mMSA/c"; do
    pred_dir=${output_dir}/${pid}/${msa}
    aln_file=${pred_dir}/${pid}.aln
    if [ -f ${pred_dir}/final.aln ]; then
      aln_file=${pred_dir}/final.aln
    fi
    if [ -f ${aln_file} ]; then
      mkdir -p ${pred_dir}/${pid}/msas
      # aln to a3m
      cat ${aln_file} | awk -v pid=${pid} 'BEGIN{i=0}{printf(">%s_%d\n", pid, i); print $0; ++i;}' > ${pred_dir}/${pid}/msas/bfd_uniref_hits.a3m
    fi
  done
else
  if [ $# -lt ${WORKER_NUM} ]; then
    WORKER_NUM=$#
  fi

  echo $@ | \
    awk '{for (i=1; i<=NF; ++i){print $i}}' | \
    ./bin/mapred -m "${deep_msa2_cmd}" -c ${WORKER_NUM}
fi
