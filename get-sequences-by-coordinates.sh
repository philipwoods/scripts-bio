#!/bin/bash

usage() {
    echo ""
    echo "get-sequences-by-coordinates.sh -t FEATURE_TABLE -f FASTA -o OUT_DIR"
    echo "  The FASTA file should contain all of the contigs referenced in the feature table."
    echo "  The feature table should not contain a header and should not have any empty lines."
    echo "  If contig names or feature names are not unique in the feature table, you may lose data."
    echo ""
}

table=""
fasta=""
output=""

while getopts ":ht:o:" opt; do
    case $opt in
        h)  usage
            exit 1
            ;;
        t)  table="${OPTARG}"
            ;;
        f)  fasta="${OPTARG}"
            ;;
        o)  output="${OPTARG}"
            ;;
        \?) echo "Invalid option: -${OPTARG}. Use the -h option for more information.">&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument. Use the -h option for more information.">&2
            exit 1
            ;;
    esac
done

if [[ ! -f "$table" ]]; then
    echo "File does not exist: $table"
    exit 1
fi
if [[ ! -f "$fasta" ]]; then
    echo "File does not exist: $fasta"
    exit 1
fi

if [[ ! -d "$output" ]]; then
    mkdir "$output"
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

temp=$(mktemp)

# Format is Contig, Feature, Start, Stop, Strand, other stuff...
# Each domain annotation has its own entry, but we only care about the overall features.
# In the first five columns, calling uniq will collapse everything to one line per feature.
cut -f 1,2,3,4,5 "$input" | uniq > $temp

while read contig name start stop strand; do
    $DIR/get-sequence-by-coordinates.py -o "$output" "$fasta" "$contig" "$name" "$start" "$stop" "$strand"
done < $temp

