#!/bin/bash

# TO DO split based on conda env requirements
# May need to make this script work on a batch instead of individual files in that case, like anvi-run-emapper

print_usage() {
    echo "bash anvi-run-fegenie.sh -c CONTIGS_DB -T NUM_THREADS -o OUT_DIR"
    echo ""
    echo ""
}

print_version() {
    echo "Last modified 14 Oct 2021"
}

contigs_db=""
out_dir=""
num_threads=""

while getopts ":hvc:o:t:" opt; do
    case $opt in
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        c)  contigs_db="${OPTARG}"
            ;;
        o)  out_dir="${OPTARG}"
            ;;
        t)  num_threads="${OPTARG}"
            ;;
        \?) echo "Invalid option: -${OPTARG}. Use the -h option for more information.">&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument. Use the -h option for more information.">&2
            exit 1
            ;;
    esac
done

if [[ ! -d ${out_dir} ]]; then
    mkdir ${out_dir}
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

filename=$(basename ${contigs_db})

echo ""
echo "  Beginning gene call export..."
echo "--------------------------------------------------------------------------------"
if [[ ! $CONDA_DEFAULT_ENV =~ "anvio" ]]; then
    echo "You need to activate an anvio conda environment before running this."
    exit 1
fi
# The parameter expansion removes the .db file extension.
gene_calls="${out_dir}/gene_calls_${filename%.*}.fa"
anvi-get-sequences-for-gene-calls -c $contigs_db -o $gene_calls --get-aa-sequences --report-extended-deflines
python $DIR/fegenie_modify_deflines.py $(realpath $gene_calls)

echo ""
echo "  Beginning FeGenie annotation..."
echo "--------------------------------------------------------------------------------"
if [[ ! $CONDA_DEFAULT_ENV =~ "fegenie" ]]; then
    echo "You need to activate the fegenie conda environment before running this."
    exit 1
fi
fegenie_out_dir="${out_dir}/fegenie_out_${filename%.*}"
fegenie_out="${fegenie_out_dir}/FeGenie-geneSummary.csv"
FeGenie.py -bin_dir "${out_dir}" -bin_ext fa -t $num_threads -out "${fegenie_out_dir}" --orfs -delim "_"

echo ""
echo "  Beginning annotation import..."
echo "--------------------------------------------------------------------------------"
if [[ ! $CONDA_DEFAULT_ENV =~ "anvio" ]]; then
    echo "You need to activate an anvio conda environment before running this."
    exit 1
fi
# For some reason, the first line with headers is doubled in this file
# With this conditional the script will keep working even if that gets fixed
if [[ $(sed -n '1p' ${fegenie_out}) == $(sed -n '2p' ${fegenie_out}) ]]; then
    sed -i 1d ${fegenie_out}
fi
anvi_table="${fegenie_out_dir}/fegenie-anvi-table.tsv"
python $DIR/fegenie_convert_table.py $(realpath $fegenie_out)

anvi-import-functions -c ${contigs_db} -i ${anvi_table}
rm -f ${gene_calls}

echo ""
echo "Done!"

