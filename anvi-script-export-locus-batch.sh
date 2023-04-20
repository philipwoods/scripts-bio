#!/bin/bash

print_usage() {
    echo "Synopsis:"
    echo "  anvi-script-export-locus-batch [-h] [-n] -i <config-input> -o <out-dir>"
    echo ""
    echo "Description:"
    echo "  This script takes in a tab-delimited configuration file and uses the information"
    echo "  to export loci from Anvi'o contigs databases in a batch. This is useful when each"
    echo "  source genome requires different export parameters."
    echo ""
    echo "  The input configuration file needs to have the following columns:"
    echo "  <genome-name>  <contigsdb>  <flank-start>  <flank-end>  <prefix>  <optional notes columns>"
    echo "  By default, the script assumes this file has a header with column names and"
    echo "  each subsequent line holds a record of a gene locus in a single genome."
    echo "  <genome-name> is a unique identifier for each genome."
    echo "  <contigsdb> is the path to the Anvi'o contigs DB for that genome."
    echo "  <flank-start> is the gene caller ID of the gene at one end of the locus."
    echo "  <flank-end> is the gene caller ID of the gene at the other end of the locus."
    echo "  <prefix> is a descriptor of the locus. Loci within a single genome cannot share descriptors."
    echo "  Further columns may be present but will not be parsed."
    echo ""
    echo "Options:"
    echo "  -h: Display this help text and exit."
    echo "  -n: Specify that the input file has no header line."
    echo "  -i: Specify the input (config file)"
    echo "  -o: Specify output directory."
}

infile=""
outdir=""
header=""

while getopts ":hni:o:" opt; do
    case $opt in 
        h)  print_usage
            exit 1
            ;;
        n)  header="noheader"
            ;;
        i)  infile="${OPTARG}"
            ;;
        o)  outdir="${OPTARG}"
            ;;
        \?) echo "Invalid option: -${OPTARG}. Use the -h option for more information.">&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument. Use the -h option for more information.">&2
            exit 1
            ;;
    esac
done

if [[ ! -e "$infile" ]]; then
    echo "File does not exist: ${infile}"
    exit 1
fi
if [[ ! -d "$outdir" ]]; then
    echo "Directory does not exist: ${outdir}"
    exit 1
fi

if [[ ! $CONDA_DEFAULT_ENV =~ "anvio" ]]; then
    echo "You need to activate an anvio conda environment before running this."
    exit 1
fi

# The input configuration file needs to have the following columns:
# <genome-name>  <contigsdb>  <flank-start>  <flank-end>  <prefix>  <optional notes columns>
# Each line holds a record of a set of ordered genes in a single genome.
# <genome-name> is a unique identifier for each genome.
# <contigsdb> is the path to the Anvi'o contigs DB for that genome.
# <flank-start> is the gene caller ID of the gene at one end of the locus.
# <flank-end> is the gene caller ID of the gene at the other end of the locus.
# <prefix> is a descriptor of the locus. Loci within a single genome cannot share descriptors.
# Any further columns may be present but will not be parsed.

data=$(mktemp)
cp $infile $data
if [[ $header == "noheader" ]]; then
    sed -i '1d' $data
fi

while read -r name db flankstart flankend locus remainder; do
    anvi-export-locus -c "$db" --flank-mode --gene-caller-ids "$flankstart","$flankend" -O "${name}_${locus}" -o "$outdir"
done < "$data"

rm -f $data

echo "Done!"

exit
