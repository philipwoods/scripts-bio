#!/bin/bash

print_usage() {
    echo "bash run-kegg-decoder.sh [-r] -i IN_DIR -o OUT_DIR"
    echo ""
    echo "  Must be run twice to complete the process. Before running for the"
    echo "  first step (extracting annotations), load the anvio-7 environment"
    echo "  in the anaconda3 module. Before running for the second step"
    echo "  (decoding the annotations), load the keggdecoder environment in"
    echo "  then anaconda3 module."
    echo ""
    echo " Step 1:"
    echo "  Takes in a directory containing Anvi'o contigs database files."
    echo "  Extracts KEGG annotations from the databases and compiles these"
    echo "  into OUT_DIR/ALL_GENES.ko. If you want to replace an existing"
    echo "  version of this file, use the -r flag to reset it."
    echo " Step 2:"
    echo "  If the ALL_GENES.ko file already exists in the output directory,"
    echo "  run KEGG-decoder on the file."
}

print_version() {
    echo "Last modified 13 Oct 2021"
}

in_dir=""
out_dir=""
reset=""

while getopts ":hvri:o:" opt; do
    case $opt in
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        r)  reset="reset"
            ;;
        i)  in_dir="${OPTARG}"
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

if [[ ! -d ${out_dir} ]]; then
    mkdir ${out_dir}
fi

all_kos="${out_dir}/ALL_GENES.ko"

# KEGG decoder takes input in the following format (no header line):
# <genome-name>_<gene-id> <KEGG-accession>
# The program uses underscores to separate genome name from gene ID, so neither of them can have underscores.
#
# anvi-export-functions outputs information in the following format:
# gene_callers_id   source  accession           ...
# <gene-id>         KOfam   <KEGG-accession>    ...
# So we have to remove the header, add the genome name to the beginning of each line, and cut just the first and third columns.

# If the output file containing compiled annotations doesn't exist,
# or we've been asked to reset it...
if [[ ! -e $all_kos || -n $reset ]]; then
    echo ""
    echo "  Extracting KEGG annotations..."
    echo "--------------------------------------------------------------------------------"
    
    # Check that the proper conda environment is loaded
    if [[ ! $CONDA_DEFAULT_ENV =~ "anvio" ]]; then
        echo "You need to load the anvio-7 conda environment before running this step."
        echo ""
        exit 1
    fi

    # Create output file, and clear any previously existing version.
    printf "" > $all_kos
    
    for i in ${in_dir}/*.db; do
        if [[ ! -e $i ]]; then
            echo "No database files in the specified directory!"
            exit 1
        fi
        filename=$(basename $i)
        # Remove extension from the file name.
        genomename=${filename%.db}
        # Replace all underscores with dashes.
        name=${genomename//_/-}
        # Pull annotations from the database
        tmp=$(mktemp)
        anvi-export-functions -c $i --annotation-sources KOfam -o $tmp
        # Remove header line from the annotation file and prepend the genome name
        sed -i '1d' $tmp # Remove header line
        sed -i "s/^/${name}_/" $tmp
        # Cut out the columns we need and append them to the output file
        cut -f 1,3 --output-delimiter=" " $tmp >> $all_kos
        rm -f $tmp
        echo "Done adding $filename"
    done
else
    echo ""
    echo "  Beginning KEGG decoding..."
    echo "--------------------------------------------------------------------------------"
    
    # Check that the proper conda environment is loaded
    if [[ ! $CONDA_DEFAULT_ENV == "keggdecoder" ]]; then
        echo "You need to load the keggdecoder conda environment before running the decoding step."
        echo "Alternatively, if you are trying to re-export KEGG annotations, you need to use the -r flag."
        echo ""
        exit 1
    fi

    decoder_out="${out_dir}/KEGG-decoder.list"
    KEGG-decoder -i $all_kos -o $decoder_out -v static
    tmp=$(mktemp)
    # Replace all dashes, slashes, or periods in the first field of the output (genome names) with underscores.
    cp "$decoder_out" "$tmp" && gawk -v FS='\t' -v OFS='\t' '{gsub(/[-\.]/,"_",$1) ; print}' $tmp > $decoder_out
    rm -f "$tmp"
    # Change the column separators to conform to a tsv format.
    sed -i 's/\t /\t/g' $decoder_out
fi
echo "Done!"

