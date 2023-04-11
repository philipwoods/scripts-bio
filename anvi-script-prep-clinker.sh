#!/bin/bash

print_usage() {
    echo "clinker-prep -i <in-dir> -o <out-dir>"
}

indir=""
outdir=""

while getopts ":hi:o:" opt; do
    case $opt in
        h)  print_usage
            exit 1
            ;;
        i)  indir="${OPTARG}"
            ;;
        o)  outdir=$(realpath ${OPTARG})
            ;;
        \?) echo "Invalid option: -${OPTARG}. Use the -h option for more information.">&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument. Use the -h option for more information.">&2
            exit 1
            ;;
    esac
done

if [[ -z "$indir" ]]; then
    echo "You must provide an input directory."
    exit 1
fi

if [[ ! -d "$indir" ]]; then
    echo "The provided directory does not exist: $indir"
    exit 1
fi

if [[ ! $CONDA_DEFAULT_ENV =~ "anvio" ]]; then
    echo "You need to activate an anvio conda environment before running this."
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

# Extract GFF of annotations
for f in ${indir}/*.db; do
    anvi-get-sequences-for-gene-calls -c $f --export-gff3 -o ${f%.db}.gff
done

# Convert GFF to GB format
for f in ${indir}/*.gff; do
    $DIR/gff_to_genbank.py $f ${f%.gff}.fa
done

# Copy GB files for clinker to the output directory
cp ${indir}/*.gb $outdir

echo "Done!"

