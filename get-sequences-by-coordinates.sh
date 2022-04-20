#!/bin/bash

usage() {
    echo ""
    echo "get-sequences-by-coordinates.sh -t FEATURE_TABLE -f FASTA -o OUT_DIR"
    echo "  The FASTA file should contain all of the contigs referenced in the feature table."
    echo "  The feature table should be formatted like a KBase domain annotation output file,"
    echo "  with each line containing the following fields in order:"
    echo ""
    echo "  Contig_name  Feature_name  Start_position  End_position  Feature_strand(+/-)"
    echo ""
    echo "  Any additional fields following these in a line will be ignored. The feature table"
    echo "  should not contain a header and should not have any empty lines. If contig names"
    echo "  or feature names are not unique in the feature table, you may lose data."
    echo "  Requires biopython to run."
    echo ""
}

table=""
fasta=""
output=""

while getopts ":ht:f:o:" opt; do
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
# Calling uniq on the first five columns will collapse everything to one line per feature.
cut -f 1,2,3,4,5 "$table" | uniq > $temp

while read contig name start stop strand; do
    $DIR/get-sequence-by-coordinates.py -o "$output" "$fasta" "$contig" "$start" "$stop" "$strand" "$name"
done < $temp

# Clean up
rm -f $temp

