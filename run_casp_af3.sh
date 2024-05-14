
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
      "  --renumber_rna_start <renumber_rna_start>\n" \
      "  --renumber_dna_start <renumber_dna_start>\n" \
      "\n"
  exit $1
}

output_dir="."
target=""
stoichiometry=""
renumber_rna_start=0
renumber_dna_start=-1

#########################
# NOTE: for MacOS
# brew install gnu-getopt
case `uname` in
    Darwin*) getopt=${getopt_path:-"/opt/homebrew/opt/gnu-getopt/bin/getopt"}
        ;;
    *) getopt=getopt
        ;;
esac

ARGS=$($getopt --options "o:h" --longoptions "output_dir:,target:,stoich:,renumber_rna_start:,renumber_dna_start:,chain_mapping:,help" -- "$@") || exit
eval "set -- ${ARGS}"
while true; do
  case "$1" in
    (-o | --output_dir) output_dir="$2"; shift 2;;
    (--target) target="$2"; shift 2;;
    (--stoich) stoichiometry="$2"; shift 2;;
    (--renumber_rna_start) renumber_rna_start="$2"; shift 2;;
    (--renumber_dna_start) renumber_dna_start="$2"; shift 2;;
    (--chain_mapping) chain_mapping="$2"; shift 2;;
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

# . ~/.local/anaconda3/bin/activate pf2

model=1
for cif in $(find $* -name "*summary_confidences_*.json" -exec grep -H ranking_score {} \; | sort -grk3 | awk -F ':' '{print $1}' | awk -F 'summary_confidences|.json' '{printf("%smodel%s.cif\n", $1, $2)}' | head -n 5); do
  echo ${cif}
  python submit_casp.py --dry_run \
    --target=${target} \
    --model=${model} \
    --stoich=${stoichiometry} \
    --renumber_chain \
    --renumber_rna_chain_start=${renumber_rna_start} \
    --renumber_dna_chain_start=${renumber_dna_start} \
    --output_file=${output_dir}/model_${model}.pdb \
    --chain_mapping="${chain_mapping}" \
    ${cif}
  model=$(expr ${model} + 1)
done
