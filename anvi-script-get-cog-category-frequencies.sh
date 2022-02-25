#/bin/bash

print_usage() {
    echo "Synopsis:"
    echo "  anvi-script-get-cog-category-frequencies [-h] [-v] -p -e EXTERNAL_GENOMES -o OUT_DIR"
    echo ""
    echo "Description:"
    echo "  To run this script, you must have already activated your anvio environment."
    echo "  Provide an Anvi'o external genomes file containing names and paths for the"
    echo "  contigs databases of interest. This script will extract the COG category"
    echo "  annotation for genes in the provided contigs databases, then compile a"
    echo "  tab-separated file of the frequency of each COG category in each genome."
    echo ""
    echo "Options:"
    echo "  -h: Display this help text and exit."
    echo "  -v: Display version text and exit."
    echo "  -p: Output in a pretty printed format."
    echo "  -e: Specify an anvio external genomes file containing genomes to use."
    echo "  -o: Specify output directory."
}

print_version() {
    echo "Last updated 19 October 2021"
}

in_file=""
out_dir=""
pretty_out="0"
pretty_suffix=""

while getopts ":hve:o:p" opt; do
    case $opt in 
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        e)  in_file="${OPTARG}"
            ;;
        o)  out_dir="${OPTARG}"
            ;;
        p)  pretty_out="1"
            pretty_suffix="-pretty"
            ;;
        \?) echo "Invalid option: -${OPTARG}. Use the -h option for more information.">&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument. Use the -h option for more information.">&2
            exit 1
            ;;
    esac
done

if [ ! -e "${in_file}" ]; then
    echo "File does not exist: ${in_file}"
    exit 1
fi
if [ ! -d "${out_dir}" ]; then
    echo "Directory does not exist: ${out_dir}"
    exit 1
fi

if [[ ! $CONDA_DEFAULT_ENV =~ "anvio" ]]; then
    echo "You need to activate an anvio conda environment before running this."
    exit 1
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Split the external genomes file into columns.
name_file=$(mktemp)
db_file=$(mktemp)
cut -f 1 "${in_file}" | sed '1d' > "${name_file}"
cut -f 2 "${in_file}" | sed '1d' > "${db_file}"

echo "Beginning COG_CATEGORY annotation export:"
# Get each (name, database) pair.
while read -u 3 genome_name && read -u 4 contigs_db; do
    # Make sure the database file exists.
    if [ ! -e "${contigs_db}" ]; then
        echo "File does not exist: ${contigs_db}"
        exit 1
    fi

    # If you change the names of any of these files, make sure to check in the
    # parse-cog-cat-functions.py helper script to make sure it will still work.
    functions_out="${out_dir}/${genome_name}_functions_COG_CATEGORY.tab"

    # Output COG_CATEGORY annotations for each database using the associated name for the output file.
    echo "Exporting annotations for ${genome_name}..."
    anvi-export-functions -c "$contigs_db" --annotation-sources COG20_CATEGORY -o "${functions_out}"
    # Cut out the 'functions' column and discard everything else.
    tmp=$(mktemp)
    cut -f 3 "${functions_out}" > "${tmp}"
    # Replace the separator '!!!' with a new line so each annotation is a single letter per line.
    sed 's/!!!/\n/g' "${tmp}" > "${functions_out}"
    # Clean up temporary files.
    rm -f "${tmp}"
done 3<"${name_file}" 4<"${db_file}"

# Helper Python script parses each file and creates a tab-delimited output.
echo "COG_CATEGORY annotation export complete."
echo "Consolidating records..."
python "${DIR}/parse-cog-cat-functions.py" "${out_dir}" $pretty_out > "${out_dir}/cog-category-frequency${pretty_suffix}.tsv"

# Clean up files
echo "Cleaning temporary files..."
while read -u 3 genome_name; do
    functions_out="${out_dir}/${genome_name}_functions_COG_CATEGORY.tab"
    rm -f "${functions_out}"
done 3<"${name_file}"

rm -f "${name_file}"
rm -f "${db_file}"

echo "Done!"

exit
