#!/bin/bash

print_usage() {
    echo "As of anvio v7.1, the built-in function anvi-compute-metabolic-enrichment does the exact thing that this script does. Use that instead."
    echo ""
    echo "anvi-script-metabolic-enrichment -m METABOLISM -g GROUPS -n COLUMN-NAME -o OUT-FILE"
    echo "  METABOLISM should be the presence-absence matrix output from anvi-estimate-metabolism"
    echo "  GROUPS should be formatted like an Anvio layers additional data file"
    echo "  COLUMN-NAME should be the name of the column in the GROUPS file that contains the group specification"
}

print_version() {
    echo "Last updated 21 October 2021"
}

metabolism_file=""
group_file=""
group_col=""
out_file=""

while getopts ":hvm:n:g:o:" opt; do
    case $opt in 
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        m)  metabolism_file="${OPTARG}"
            ;;
        n)  group_col="${OPTARG}"
            ;;
        g)  group_file="${OPTARG}"
            ;;
        o)  out_file="${OPTARG}"
            ;;
        \?) echo "Invalid option: -${OPTARG}. Use the -h option for more information.">&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument. Use the -h option for more information.">&2
            exit 1
            ;;
    esac
done

if [[ ! -e $metabolism_file ]]; then
    echo "File does not exist: $metabolism_file"
    exit 1
fi
if [[ ! -e $group_file ]]; then
    echo "File does not exist: $group_file"
    exit 1
fi
if [[ -z $group_col ]]; then
    echo "You must provide a column name in the group file to pull the groups from."
    exit 1
fi

if [[ ! $CONDA_DEFAULT_ENV =~ 'anvio' ]]; then
    echo "You need to activate an anvio conda environment before running this."
    exit 1
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo "Formatting the metabolism matrix file..."
tmp=$(mktemp)
python "${DIR}/format-enrichment-input.py" -m "${metabolism_file}" -n "${group_col}" -g "${group_file}" > $tmp
echo "Computing enrichment..."
$CONDA_PREFIX/bin/anvi-script-enrichment-stats --input=$tmp --output=$out_file
echo "Done!"
