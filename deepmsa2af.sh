#!/bin/bash
#

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

help() {
  echo -e "usage: `basename $0` [options] input_fasta\n" \
          "options:\n" \
          "  -o, --deepmsa2_output_dir <deepmsa2_output_dir>\n" \
          "      --hhr_file <hhr_file>\n"
  exit $1
}

deepmsa2_output_dir="."
hhr_file=""

ARGS=$(getopt --options "o:h" --longoptions "deepmsa2_output_dir:,hhr_file:,help" -- "$@") || exit
eval "set -- ${ARGS}"
echo "ARGS: $*"
while true; do
  case "$1" in
    (-o | --deepmsa2_output_dir) deepmsa2_output_dir="$2"; shift 2;;
    (--hhr_file) hhr_file="$2"; shift 2;;
    (-h | --help) help 0 ;;
    (--) shift 1; break;;
    (*) help 1;
  esac
done
echo "deepmsa2_output_dir: ${deepmsa2_output_dir}"

input_fasta=$1;
if [ -z "${input_fasta}" ]; then
  help 1;
fi
if [ -z "${hhr_file}" ]; then
  echo "WARNING: You are makeing AF2 format MSAs without template file !!!"
fi

fasta_name=$(basename ${input_fasta}|awk -F '.' '{printf("%s", $1); for(i=2;i<NF;++i) printf(".%s", $i)}')

for c in "dMSA" "qMSA" "mMSA/a" "mMSA/b" "mMSA/c"; do
  output_dir=${deepmsa2_output_dir}/${c}/${fasta_name}
  aln=${output_dir}.aln
  if [ -e ${aln} ]; then
    mkdir -p ${output_dir}/msas
    cat ${aln} | awk -v n=${fasta_name} 'BEGIN{i=0}{printf(">%s_%d\n", n, i); print $0; ++i;}' > ${output_dir}/msas/bfd_uniref_hits.a3m
    if [ ! -z "${hhr_file}" ]; then
      cp ${hhr_file} ${output_dir}/msas/pdb_hits.hhr
    fi
  else
    echo "WARNING: ${aln} don't exists!"
  fi
done
