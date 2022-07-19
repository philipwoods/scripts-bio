#!/bin/bash

print_usage() {
    echo ""
    echo "anvi-script-frequency-enrichment -e <external-genomes> -g <group-file> -a <annotation> -m <mode> -o <out-dir>"
    echo "  Performs an enrichment analysis comparable to anvi-metabolic-enrichment or anvi-compute-functional-enrichment"
    echo "  except it analyzes the frequency of each function in each genome of a group rather than simply its presence"
    echo "  or absence in each genome of a group."
    echo ""
    echo "  <external-genomes> should be formatted like an Anvio external genomes file."
    echo "  <group-file> should be formatted like an Anvio layers additional data file with a column labeled 'group' specifying group membership."
    echo "  <annotation> is the annotation source to use when analyzing frequency enrichment."
    echo "  <mode> is the output mode to use when determining annotation frequency. The available modes are specified below:"
    echo "      annotation  :   The frequency is the proportion of annotations in the group with a given value."
    echo "      gene        :   The frequency is the proportion of genes in the group with a given annotation."
    echo ""
}

print_version() {
    echo "Last updated 4 November 2021"
}

in_file=""
group_file=""
out_dir=""
annotation=""
mode=""

while getopts ":hve:g:a:m:o:" opt; do
    case $opt in 
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        e)  in_file="${OPTARG}"
            ;;
        g)  group_file="${OPTARG}"
            ;;
        a)  annotation="${OPTARG}"
            ;;
        m)  mode="${OPTARG}"
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
if [[ -z $mode ]]; then
    echo "You must select an output mode. Use the -h option for more information." >&2
    exit 1
fi

if [[ ! $CONDA_DEFAULT_ENV =~ 'anvio' ]]; then
    echo "You need to activate an anvio conda environment before running this." >&2
    exit 1
fi

if [[ -z "${annotation}" ]]; then
    annotation="COG20_CATEGORY"
    echo "WARNING: You chose not to provide an annotation source, so ${annotation} will be used by default."
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

# Split the external genomes file into columns.
name_file=$(mktemp)
db_file=$(mktemp)
cut -f 1 "${in_file}" | sed '1d' > "${name_file}"
cut -f 2 "${in_file}" | sed '1d' > "${db_file}"

echo "Beginning annotation export:"
# Get each (name, database) pair.
while read -u 3 genome_name && read -u 4 contigs_db; do
    # Make sure the database file exists.
    if [ ! -e "${contigs_db}" ]; then
        echo "File does not exist: ${contigs_db}"
        exit 1
    fi

    # Output annotations for each database using the associated name for the output file.
    echo "Exporting annotations for ${genome_name}..."
    # As of anvio v7.1, the function export file is formatted as follows:
    # gene_callers_id   source   accession   function   e_value
    tmp=$(mktemp)
    anvi-export-functions -c "$contigs_db" --annotation-sources "${annotation}" -o "${tmp}"
    # Cut out the 'function' column and discard everything else.
    # If you change the name of the <functions-out> file, make sure to check in the
    # parse-function-frequency.py helper script to make sure it will still work.
    functions_out="${out_dir}/${genome_name}-functions.tmp"
    cut -f 3,4  "${tmp}" > "${functions_out}"
    # Clean up temporary files.
    rm -f "${tmp}"
done 3<"${name_file}" 4<"${db_file}"

# Helper Python script parses each file and creates tab-delimited output.
echo "Annotation export complete."
echo "Consolidating records..."
freq_file="${out_dir}/frequency-${annotation}-${mode}.tsv"
counts_file="${out_dir}/${mode}-counts.tmp"
python "${DIR}/parse-function-frequency.py" "${out_dir}" "${mode}" "${freq_file}" "${counts_file}"
# The counts file created by the helper script is useful in annotation mode,
# but in gene mode we need to make some changes.
if [[ $mode == "gene" ]]; then
    echo "Counting genes..."
    printf "genome\tN\n" > "${counts_file}"
    while read -u 3 genome_name && read -u 4 contigs_db; do
        tmp=$(mktemp)
        anvi-export-gene-calls -c "${contigs_db}" -o $tmp --gene-caller prodigal --skip-sequence-reporting
        # The gene calls file has one line per gene call and one additional line for the header.
        count=$(wc -l < $tmp)
        count=$((count-1))
        printf "${genome_name}\t${count}\n" >> "${counts_file}"
        rm -f $tmp
    done 3<"${name_file}" 4<"${db_file}"
fi

# Clean up files
echo "Cleaning temporary files..."
while read -u 3 genome_name; do
    functions_out="${out_dir}/${genome_name}-functions.tmp"
    rm -f "${functions_out}"
done 3<"${name_file}"

rm -f "${name_file}"
rm -f "${db_file}"

###

echo "Formatting the function frequency file..."
tmp=$(mktemp)
python "${DIR}/format-enrichment-input.py" -f "${freq_file}" -g "${group_file}" -c "${counts_file}" > $tmp
rm -f "${counts_file}"
echo "Computing enrichment..."
out_file="${out_dir}/enrichment-${annotation}-${mode}.tsv"
$CONDA_PREFIX/bin/anvi-script-enrichment-stats --input=$tmp --output=$out_file

echo "Done!"

