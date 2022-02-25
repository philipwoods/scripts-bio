#!/bin/bash

print_usage() {
    echo "Synopsis:"
    echo "  anvi-script-get-sequences-by-function [-h|-v] [-a] -f ANNOTATION_SOURCE"
    echo "                                        -s SEARCH -e GENOMES_FILE -o OUT_DIR"
    echo ""
    echo "Description:"
    echo "  To run this script, you must have activated your anvio environment."
    echo "  This script takes in an external genomes file as used by Anvi'o and searches"
    echo "  the associated genomes for genes based on functional annotation. These gene"
    echo "  sequences are then exported into a FASTA file either as nucleotide or amino"
    echo "  acid sequences, depending on the use of the -a flag. Nucleotide sequences are"
    echo "  exported by default."
    echo ""
    echo "  Assumes at most one annotation per gene will match the search term. Because"
    echo "  of this, the EGGNOG_BACT, COG20_FUNCTION, and KOfams annotation sources are"
    echo "  well suited to this analysis. Using Pfams may be possible, but may fail."
    echo ""
    echo "Options:"
    echo "  -h: Display this help text and exit."
    echo "  -v: Display version text and exit."
    echo "  -a: Export amino acid sequences instead of nucleotide sequences."
    echo "  -e: Specify Anvi'o external genomes file to use."
    echo "  -f: Specify annotation source for the function search."
    echo "  -s: Specify regex search term to use for selecting functions."
    echo "  -o: Specify output directory."
}

print_version() {
    echo "Last modified 20 February 2020"
}

genome_file=""
search_term=""
annotation_src=""
out_dir=""
aa_option=""

while getopts ":hve:s:f:o:a" opt; do
    case $opt in
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        e)  genome_file="${OPTARG}"
            ;;
        s)  search_term="${OPTARG}"
            ;;
        f)  annotation_src="${OPTARG}"
            ;;
        o)  out_dir="${OPTARG}"
            ;;
        a)  aa_option="--get-aa-sequences"
            ;;
        \?) echo "Invalid option: -${OPTARG}. Use the -h option for more information.">&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument. Use the -h option for more information.">&2
            exit 1
            ;;
    esac
done

if [[ ! -e "${genome_file}" ]]; then
    echo "File does not exist: ${genome_file}"
    exit 1
fi
if [[ ! -d "${out_dir}" ]]; then
    echo "Directory does not exist: ${out_dir}"
    exit 1
fi
if [[ ! $CONDA_DEFAULT_ENV =~ "anvio" ]]; then
    echo "You need to activate an anvio conda environment before running this."
    exit 1
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

functions_temp="${out_dir}/functions_temp.tab"
sequences_temp="${out_dir}/sequence_temp.fasta"

sequences_out="${out_dir}/${annotation_src}_${search_term}.fasta"
printf "" > "${sequences_out}"

while read name db; do
    echo "Searching for ${search_term} in ${name}..."
    anvi-export-functions -c "${db}" --annotation-sources "${annotation_src}" -o "${functions_temp}"
    gene_id=$(grep "${search_term}" "${functions_temp}" | cut -f 1 | tr '\n' ',' | sed 's/.$//')
    gene_info=$(grep "${search_term}" "${functions_temp}" | cut -f 2 --complement | sed 's/\t/|/g')
    if [[ -z "${gene_id}" ]]; then
        echo "No results found."
    else
        echo "Exporting gene sequences from ${name}..."
        if [[ -z "${aa_option}" ]]; then
            anvi-get-sequences-for-gene-calls -c "${db}" --gene-caller-ids "${gene_id}" -o "${sequences_temp}"
        elif [[ -n "${aa_option}" ]]; then
            anvi-get-sequences-for-gene-calls -c "${db}" --gene-caller-ids "${gene_id}" -o "${sequences_temp}" "${aa_option}"
        fi
        python "${DIR}/modify_fasta_defline.py" "${sequences_temp}" "${name}" <(echo "${gene_info}")
        cat "${sequences_temp}" >> "${sequences_out}"
        echo "" >> "${sequences_out}"
        rm -f "${sequences_temp}"
    fi 
    rm -f "${functions_temp}"
done < <(sed 1d "${genome_file}")
echo "Done!"
exit
