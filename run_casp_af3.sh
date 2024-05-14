
#!/bin/bash
#

set -e

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

cd ${CWD}

help() {
  echo -e "usage: `basename $0` [options] af3_fold_dir1 af3_fold_dir2 ...\n" \
          "options:\n" \
	  "  -o, --output_dir <output_dir> *required*\n" \
	  "  --target <target> *required* ex: T1201\n" \
	  "  --stoich <stoichiometry> *required*\n" \
	  "\n"
  exit $1
}

output_dir="."
target=""
stoichiometry=""

ARGS=$(getopt --options "o:h" --longoptions "output_dir:,target:,stoich:,help" -- "$@") || exit
eval "set -- ${ARGS}"
while true; do
  case "$1" in
    (-o | --output_dir) output_dir="$2"; shift 2;;
    (--target) target="$2"; shift 2;;
    (--stoich) stoichiometry="$2"; shift 2;;
    (-h | --help) help 0 ;;
    (--) shift 1; break;;
    (*) help 1;
  esac
done

if [ -z "${stoichiometry}" ]; then
  help 1;
fi
if [ -z "${target}" ]; then
  help 1;
fi

mkdir -p ${output_dir}

model=1
for cif in $(find $* -name *summary_confidences_*.json -exec grep -H ranking_score {} \; | sort -grk3 | awk -F ':' '{print $1}' | awk -F 'summary_confidences|.json' '{printf("%smodel%s.cif\n", $1, $2)}' | head -n 5); do
  echo ${cif}
  python submit_casp.py --dry_run -v \
    --target=${target} \
    --model=${model} \
    --stoich=${stoichiometry} \
    --renumber_chain \
    --output_file=${output_dir}/model_${model}.pdb \
    ${cif}
  model=$(expr ${model} + 1)
done
