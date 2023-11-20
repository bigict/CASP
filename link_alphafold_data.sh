#!/bin/bash
#

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

cd ${CWD}

help() {
  echo "usage: `basename $0` -s <source_data_dir> -t <target_data_dir> -f <pdb70_unlink_file>"
  exit $1
}

pdb70_unlink_file=""
while getopts 's:t:f:h' OPT; do
  case $OPT in
    s) source_data_dir="$OPTARG";;
    t) target_data_dir="$OPTARG";;
    f) pdb70_unlink_file="$OPTARG";;
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

for d in "params" "pdb70"; do
  echo "cp -r ${source_data_dir}/${d} ${target_data_dir}/${d}"
  cp -r ${source_data_dir}/${d} ${target_data_dir}/${d}
done

for d in "params" "uniref90" "mgnify" "bfd" "uniclust30" "pdb_mmcif"; do
  echo "ln -s -f ${source_data_dir}/${d} ${target_data_dir}/${d}"
  ln -s -f ${source_data_dir}/${d} ${target_data_dir}/${d}
done

db="pdb70"
if [ ! -z ${pdb70_unlink_file} ]; then
  chmod u+w ${target_data_dir}/${db}/${db}_a3m.ffindex
  ffindex_modify -s -u -f ${pdb70_unlink_file} ${target_data_dir}/${db}/${db}_a3m.ffindex
  chmod u+w ${target_data_dir}/${db}/${db}_hhm.ffindex
  ffindex_modify -s -u -f ${pdb70_unlink_file} ${target_data_dir}/${db}/${db}_hhm.ffindex
  chmod u+w ${target_data_dir}/${db}/${db}_cs219.ffindex
  ffindex_modify -s -u -f ${pdb70_unlink_file} ${target_data_dir}/${db}/${db}_cs219.ffindex
fi
