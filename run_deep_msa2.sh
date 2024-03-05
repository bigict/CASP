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

if [ $# -eq 0 ]; then
  while read input_fasta; do
    bash $0 -x ${CONDA_ENV} -n 1 -o ${output_dir} ${input_fasta}
  done
elif [ $# -eq 1 ]; then
  pid=$(basename $1 | awk -F '.' '{printf("%s", $1); for(i=2;i<NF;++i){printf(".%s", $i);}}');
  bash -x msa-suite/examples/deepMSA2.sh -o ${output_dir}/${pid} -i $1 -vk
else
  if [ $# -lt ${WORKER_NUM} ]; then
    WORKER_NUM=$#
  fi

  echo $@ | \
    awk '{for (i=1; i<=NF; ++i){print $i}}' | \
    ./bin/mapred -m "${deep_msa2_cmd}" -c ${WORKER_NUM}
fi
