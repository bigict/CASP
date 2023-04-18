#!/bin/bash
#

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

cd ${CWD}

help() {
  echo "usage: `basename $0` -s <source_data_dir> -t <target_data_dir>"
  exit $1
}

while getopts 's:t:h' OPT; do
  case $OPT in
    s) source_data_dir="$OPTARG";;
    t) target_data_dir="$OPTARG";;
    h) help 0;;
    ?) help 1;;
  esac
done
shift $((OPTIND - 1))

if [ -z ${source_data_dir} ]; then
  exit 1;
fi
if [ -z ${target_data_dir} ]; then
  exit 1;
fi

source_data_dir=$(realpath ${source_data_dir})
target_data_dir=$(realpath ${target_data_dir})

echo ${source_data_dir}
echo ${target_data_dir}

mkdir -p ${target_data_dir}

for d in "params" "uniprot" "uniref90" "mgnify" "bfd" "uniref30" "pdb70" "pdb_mmcif" "pdb_seqres" "pdb_seqres.txt"; do
  echo ${d}
  ln -s -f ${source_data_dir}/${d} ${target_data_dir}/${d}
done
