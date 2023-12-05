#!/bin/bash
#

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

cd ${CWD}

help() {
  echo "usage: `basename $0` -d <data_dir> -o <output_dir> -p <msa_params> input_fasta ..."
  exit $1
}

msa_tool_home="/public/dbu"
msa_tool_srv="dbu@10.206.0.2"
msa_tool_sync="rsync -avP -e ssh"
msa_tool_run="ssh -i /home/tanqian/.ssh/id_rsa ${msa_tool_srv}"

${msa_tool_run} 

#timestamp=$(date +%Y%m%d-%H%M%S-%N)
timestamp= "test"

msa_tool_data_dir="${msa_tool_home}/CASP/data/${timestamp}"
${msa_tool_run} mkdir -p ${msa_tool_data_dir}/fasta
#fasta上传远程
${msa_tool_sync} $* ${msa_tool_srv}:${msa_tool_data_dir}/fasta/

#远程对input1_fasta input2_fasta ...run SeqParser.py
for f in $@; do
  f=$(basename ${f})
  ${msa_tool_run} python ${msa_tool_home}/CASP/FUpred/SeqParse.py -i ${msa_tool_data_dir}/fasta/${f} -o ${msa_tool_data_dir}/fasta
done

echo "parse_domain done"

exit $1

af2_pred_data_dir="${CWD}/data/${timestamp}"
mkdir -p ${af2_pred_data_dir}

#将远程的${msa_tool_data_dir}/fasta分好的同步到本地
parse_data_dir="${af2_pred_data_dir}/parse"
mkdir -p ${parse_data_dir}
${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/fasta/* ${parse_data_dir}/


#创建msa_tool_file_list文件，里面存放的是各个fasta文件,空格分隔
msa_tool_file_list=""
for f in "${parse_data_dir}"/fasta/* ; do
  f=$(basename ${f})
  msa_tool_file_list="${msa_tool_file_list} ${msa_tool_data_dir}/fasta/${f}"
done

${msa_tool_run} bash -x ${msa_tool_home}/CASP/run_msa_tool.sh \
    -d ${msa_tool_home}/data/af2_db \
    -o ${msa_tool_data_dir}/msa \
    ${msa_tool_file_list}

echo "msa_tool_run done"

#check whether exists msa ,use msa_tool_file_list
for f in $@;do
  f=$(basename ${f})
  f=${f%.*}
  ${msa_tool_run} bash -x ${msa_tool_home}/CASP/Splice.sh \
    -i "${msa_tool_data_dir}/msa" \
    -n "${f}" \
    -o "${msa_tool_data_dir}/msa/${f}_splice"
done

exit $1
${msa_tool_sync} ${msa_tool_srv}:${msa_tool_data_dir}/msa/* ${af2_pred_data_dir}/

# bash -x run_alphafold.sh \
#     -m alphafold \
#     -d db \
#     -o ${af2_pred_data_dir} \
#     -p "--use_precomputed_msas" \
#     $*
