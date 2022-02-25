#!/bin/bash

print_usage() {
    echo ""
    echo "Usage"
    echo "  add-info-to-fasta.sh [-h] [-v] -c CONTIGSDB -f FASTA -a ANNOTATION_SOURCE"
    echo ""
    echo "Only one annotation source is allowed at a time. The CONTIGSDB is from"
    echo "Anvi'o, and the FASTA should be a file created from the same CONTIGSDB"
    echo "using e.g. anvi-get-sequences-for-gene-calls."
    echo ""
}

print_version() {
    echo "Last modified 18 September 2020"
}

fasta=""
contigsdb=""
annotation=""

while getopts ":hvc:f:a:" opt; do
    case $opt in
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        c)  contigsdb=$(realpath ${OPTARG})
            ;;
        f)  fasta=$(realpath ${OPTARG})
            ;;
        a)  annotation="${OPTARG}"
            ;;
        \?) echo "Invalid option: -${OPTARG}. Use the -h option for more information.">&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument. Use the -h option for more information.">&2
            exit 1
            ;;
    esac
done

if [[ -z $annotation ]]; then
    echo "You must provide an annotation source to query your contigs database."
    exit 1
fi
if [[ $annotation =~ "," ]]; then
    echo "Only one annotation source is allowed at a time."
    exit 1
fi
if [[ ! -e $contigsdb ]]; then
    echo "The provided file does not exist: $contigsdb"
    exit 1
fi
if [[ ! -e $fasta ]]; then
    echo "The provided file does not exist: $fasta"
    exit 1
fi

if [[ ! $CONDA_DEFAULT_ENV =~ "anvio" ]]; then
    echo "You need to load an anvio conda environment before running this script."
    exit 1
fi

funcs=$(mktemp)
anvi-export-functions -c $contigsdb --annotation-sources $annotation -o $funcs
tmp=$(mktemp)
cut -f 2 --complement $funcs | sed 's/\t/|/g' | sed '1d' > $tmp
mv $tmp $funcs
rm -f $tmp
python /data1/projects/pwoods/scripts/modify_fasta_defline.py $fasta $(basename $contigsdb .db) $funcs
rm -f $funcs

