#!/bin/bash

print_usage() {
    echo "Synopsis:"
    echo "  anvi-script-get-sequences-by-id [-h|-v] [-a] -i ID_FILE -e GENOMES_FILE -o OUT_DIR"
    echo ""
    echo "Description:"
    echo "  To run this script, you must have activated your anvio environment."
    echo "  This script takes in an external genomes file as used by Anvi'o and searches"
    echo "  the associated genomes for genes based on gene caller IDs. These gene"
    echo "  sequences are then exported into a FASTA file either as nucleotide or amino"
    echo "  acid sequences, depending on the use of the -a flag. Nucleotide sequences are"
    echo "  exported by default. The ID file should have two tab-separated columns: id"
    echo "  and function. The id column should have a comma-separated list of gene IDs"
    echo "  (i.e. 1,2,3), and the function column should have a comma-separated list of"
    echo "  notes corresponding to the ids (e.g. Function A,Function B,Function C). Each"
    echo "  gene ID must have a corresponding note. If you don't want any genes from a"
    echo "  particular genome, leave that line of the ID file blank."
    echo ""
    echo "Options:"
    echo "  -h: Display this help text and exit."
    echo "  -v: Display version text and exit."
    echo "  -a: Export amino acid sequences instead of nucleotide sequences."
    echo "  -e: Specify Anvi'o external genomes file to use."
    echo "  -i: Specify file listing gene IDs corresponding to each genome."
    echo "  -o: Specify output directory."
}

print_version() {
    echo "Last modified 20 February 2020"
}

genome_file=""
id_file=""
out_dir=""
aa_option=""

while getopts ":hve:i:o:a" opt; do
    case $opt in
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        e)  genome_file="${OPTARG}"
            ;;
        i)  id_file="${OPTARG}"
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
if [[ ! -e "${id_file}" ]]; then
    echo "File does not exist: ${id_file}"
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

# Get real directory where script file is located,
# even if the script was called via a symlink
SOURCE=${BASH_SOURCE[0]}
while [ -L "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )
  SOURCE=$(readlink "$SOURCE")
  [[ $SOURCE != /* ]] && SOURCE=$DIR/$SOURCE # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )

sequences_temp="${out_dir}/sequence_temp.fasta"

id_file_base=$(basename "${id_file}" ".txt")
sequences_out="${out_dir}/${id_file_base}_sequences.fasta"
printf "" > "${sequences_out}"

while read name db ids functions; do
    if [[ -z "${ids}" ]]; then
        echo "No genes requested from ${name}. Moving on..."
    else
        echo "Exporting gene sequences from ${name}..."
        if [[ -z "${aa_option}" ]]; then
            anvi-get-sequences-for-gene-calls -c "${db}" --gene-caller-ids "${ids}" -o "${sequences_temp}"
        elif [[ -n "${aa_option}" ]]; then
            anvi-get-sequences-for-gene-calls -c "${db}" --gene-caller-ids "${ids}" -o "${sequences_temp}" "${aa_option}"
        fi
        gene_info="$( paste -d '|' <(echo ${ids} | sed 's/,/\n/g') <(echo ${functions} | sed 's/,/\n/g') )"
        python $DIR/modify_fasta_defline.py "${sequences_temp}" "${name}" <(echo "${gene_info}")
        cat "${sequences_temp}" >> "${sequences_out}"
        echo "" >> "${sequences_out}"
        rm -f "${sequences_temp}"
    fi
done < <( paste "${genome_file}" "${id_file}" | sed 1d )
echo "Done!"
exit
