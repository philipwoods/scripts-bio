#!/bin/bash

print_usage() {
    echo ""
    echo "anvi-script-alignment-enrichment -f <fasta-alignment> -g <group-file> -o <out-dir>"
    echo "  Performs an enrichment analysis comparable to anvi-metabolic-enrichment or anvi-compute-functional-enrichment"
    echo "  except it analyzes the frequency of amino acids in each column of an alignment by group."
    echo ""
    echo "  <fasta-alignment> should be an aligned protein fasta file."
    echo "  <group-file> should be formatted like an Anvio layers additional data file with a column labeled 'group' specifying group membership."
    echo ""
}

print_version() {
    echo "Last updated 6 June 2024"
}

in_file=""
group_file=""
out_dir=""
annotation=""
mode=""

while getopts ":hvf:g:o:" opt; do
    case $opt in 
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        f)  in_file="${OPTARG}"
            ;;
        g)  group_file="${OPTARG}"
            ;;
        o)  out_dir="${OPTARG}"
            ;;
        \?) echo "Invalid option: -${OPTARG}. Use the -h option for more information.">&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument. Use the -h option for more information.">&2
            exit 1
            ;;
    esac
done

if [[ ! -e $in_file ]]; then
    echo "File does not exist: $in_file" >&2
    exit 1
fi
if [[ ! -e $group_file ]]; then
    echo "File does not exist: $group_file" >&2
    exit 1
fi

if [[ ! $CONDA_DEFAULT_ENV =~ 'anvio' ]]; then
    echo "You need to activate an anvio conda environment before running this." >&2
    exit 1
fi

# Get real directory where script file is located,
# even if the script was called via a symlink
SOURCE=${BASH_SOURCE[0]}
while [ -L "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )
  SOURCE=$(readlink "$SOURCE")
  [[ $SOURCE != /* ]] && SOURCE=$DIR/$SOURCE # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )
###

## PSEUDOCODE
# rewrite parse-function-frequency into parse-alignment-column-frequencies
# counts are alignment lengths or protein length? mode switch?


# Helper Python script parses each file and creates tab-delimited output.
freq_file="${out_dir}/residue-frequencies.tsv"
counts_file="${out_dir}/residue-counts.tmp"
python "${DIR}/parse-alignment-column-frequencies.py" "${out_dir}" "${freq_file}" "${counts_file}"

###

echo "Formatting the function frequency file..."
tmp=$(mktemp)
python "${DIR}/format-enrichment-input.py" -f "${freq_file}" -g "${group_file}" -c "${counts_file}" > $tmp
rm -f "${counts_file}"
echo "Computing enrichment..."
out_file="${out_dir}/enrichment-${annotation}-${mode}.tsv"
$CONDA_PREFIX/bin/anvi-script-enrichment-stats --input=$tmp --output=$out_file

echo "Done!"

