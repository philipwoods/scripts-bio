#!/bin/bash

print_usage() {
    echo "anvi-script-extract-prot-motifs [-h] [-v] -e <ext-genomes> -s <search-motif> -o <out-dir>"
    echo ""
    echo "  Searches for the provided protein motif in the genomes listed in the external genomes file."
    echo "  This file must be formatted as an Anvi'o external genomes file. The search motif should be"
    echo "  formatted as a regular expression e.g. as C..CH instead of CXXCH."
    echo ""
    echo "  <ext-genomes>  : An Anvi'o external genomes file listing genome names and contigs databases."
    echo "  <search-motif> : The protein motif to look for. Format as a regular expression."
    echo "  <out-dir>      : The directory to write output files into."
}

print_version() {
    echo "Last modified 25 November 2021"
}

summarize_counts() {
# Usage: summarize_counts <countfile> <summary>
# <countfile> is the output file produced by prot_motif_count.pl
# <summary> is a path for an output file
#
# <summary> will be formatted with two columns. The first column will contain all of the
# unique values for how many times the motif was detected within a single protein. The
# second column will contain how many unique proteins had that number of motifs present.
    local count=$1
    local summary=$2
    echo -e "motif_count\tfrequency" > "$summary"
    sed 1d "$count" | cut -f 2 | sort -g | uniq -c | awk 'OFS="\t" {print $2, $1}' >> "$summary"
}

ext_genomes=""
motif=""
outdir=""

while getopts ":hve:s:o:" opt; do
    case $opt in
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        e)  ext_genomes="${OPTARG}"
            ;;
        s)  motif="${OPTARG}"
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

if [[ ! -e "$ext_genomes" ]]; then
    echo "File does not exist: ${ext_genomes}"
    exit 1
fi
if [[ -z "$motif" ]]; then
    echo "Please provide a protein motif to search for."
    exit 1
fi

if [[ ! $CONDA_DEFAULT_ENV =~ "anvio" ]]; then
    echo "You need to activate an anvio conda environment before running this."
    exit 1
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Create output directories if necessary.
if [[ ! -d "${outdir}" ]]; then
    mkdir "${outdir}"
fi
countsdir="${outdir}/counts"
searchdir="${outdir}/search"
if [[ ! -d "${countsdir}" ]]; then
    mkdir "${countsdir}"
fi
if [[ ! -d "${searchdir}" ]]; then
    mkdir "${searchdir}"
fi

# Split the external genomes file into columns and remove the header.
name_file=$(mktemp)
db_file=$(mktemp)
cut -f 1 "${ext_genomes}" | sed '1d' > "${name_file}"
cut -f 2 "${ext_genomes}" | sed '1d' > "${db_file}"

# Get each (name, database) pair.
while read -u 3 genome_name && read -u 4 contigs_db; do
    # Make sure the database file exists.
    if [ ! -e "${contigs_db}" ]; then
        echo "File does not exist: ${contigs_db}"
        exit 1
    fi

    # Export the protein sequences for the genome and search for the motif.
    echo "Getting protein sequences for ${genome_name}..."
    proteins=$(mktemp)
    anvi-get-sequences-for-gene-calls -c "$contigs_db" --get-aa-sequences --report-extended-deflines -o "$proteins"
    echo "Searching for ${motif} in ${genome_name} proteins..."
    # The following perl functions append to their output files without
    # overwriting so we need to clear them in case they already exist.
    countsfile="${countsdir}/${genome_name}.tsv"
    searchfile="${searchdir}/${genome_name}.tsv"
    countsummary="${countsfile%.tsv}-summary.tsv"
    printf "" > ${countsfile}
    printf "" > ${searchfile}
    $DIR/prot_motif_count.pl "$proteins" "$motif" "${countsfile}"
    $DIR/prot_motif_search.pl "$proteins" "$motif" "${searchfile}"
    summarize_counts "${countsfile}" "${countsummary}"

    # Clean up temporary files.
    rm -f "${proteins}"
done 3<"${name_file}" 4<"${db_file}"
echo "Search completed."

# Create overall summary file using python for simplicity.
echo "Creating summary files..."
python $DIR/consolidate-series.py --tail="-summary.tsv" --fill-index --fill-data 0 --index "Motif count"  ${countsdir}/*-summary.tsv > "${outdir}/summary-frequency.tsv"
python $DIR/consolidate-series.py --tail="-summary.tsv" --fill-index --fill-data 0 --index "Motif count" -c ${countsdir}/*-summary.tsv > "${outdir}/summary-cumulative.tsv"
# Create info file for reference
info="${outdir}/info.txt"
printf "" > $info
echo "Search motif: ${motif}" >> $info
echo "Genomes file: $(realpath ${ext_genomes})" >> $info

echo "Done!"

