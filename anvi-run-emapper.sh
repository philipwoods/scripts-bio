#!/bin/bash

print_usage() {
    echo ""
    echo "bash anvi-run-emapper.sh -i IN_DIR -r RESULTS_DIR [-e EMAPPER_VERSION]"
    echo "  The input directory should contain anvio contigs databases to annotate."
    echo "  Other files will be ignored. Use the -e option after running emapper."
    echo ""
}

print_version() {
    echo "Last modified 25 Jul 2024"
}

in_dir=""
results_dir=""
emapper_version=""

while getopts ":hvi:r:e:" opt; do
    case $opt in
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        i)  in_dir="${OPTARG}"
            ;;
        r)  results_dir="${OPTARG}"
            ;;
        e)  emapper_version="${OPTARG}"
            ;;
        \?) echo "Invalid option: -${OPTARG}. Use the -h option for more information.">&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument. Use the -h option for more information.">&2
            exit 1
            ;;
    esac
done

if [[ ! $CONDA_DEFAULT_ENV =~ "anvio" ]]; then
    echo "You need to load an anvio conda environment before running this."
    exit 1
fi

if [[ -z $emapper_version ]]; then # If we are running this before emapper...
    if [[ ! -d ${results_dir} ]]; then
        mkdir ${results_dir}
    fi

    echo ""
    echo "  Beginning gene call export and concatenation..."
    echo "--------------------------------------------------------------------------------"

    all_genes="${results_dir}/ALL_GENES.fa"

    for i in ${in_dir}/*.db; do
        if [[ ! -e $i ]]; then
            echo "No contigs database files in the specified input directory!"
            exit 1
        fi
        filename=$(basename $i)
        # The parameter expansion removes the .db file extension.
        gene_calls="${results_dir}/gene_calls_${filename%.*}.fa"
        anvi-get-sequences-for-gene-calls -c $i -o $gene_calls --get-aa-sequences
        sed "s/^>/>${filename}|/" $gene_calls >> $all_genes
        rm -f $gene_calls
        echo "Done adding $filename"
    done

    echo ""
    echo "  NEXT STEPS"
    echo "--------------------------------------------------------------------------------"
    echo "Load the eggnog_mapper environment in anaconda3."
    echo "Verify the file ALL_GENES.fa is in your results directory, then run the following:"
    echo ""
    echo "emapper.py --data_dir /data1/db/emapper/ -m diamond -o ALL_GENES -i ${results_dir}/ALL_GENES.fa --output_dir ${results_dir} --cpu <num_threads>"
    echo ""
    echo "Once this completes, run emapper.py --version and note the version number."
    echo "Rerun this script, pointing to the same directories as the first time, and also providing"
    echo "the -e option with the appropriate version number."
else # If we are running this after emapper...
    results_file="${results_dir}/ALL_GENES.emapper.annotations"
    if [[ ! -e $results_file ]]; then
        echo "EggNOG mapper annotations file is missing from the results directory."
        echo "Are you sure you already ran emapper?"
        exit 1
    fi
    
    echo ""
    echo "  Splitting annotation results by genome..."
    echo "--------------------------------------------------------------------------------"
    tmp=$(mktemp)
    grep '#' $results_file > "$tmp"
    # Clear any old files present to avoid bugs with the >> below.
    for i in for i in $results_dir/*.emapper.annotations; do
        if [[ $(basename $i) != "ALL_GENES.emapper.annotations" ]]; then 
            rm -f $i
        fi
    done
    # Use the genome name of the query to split it into the appropriate file.
    while IFS="|" read -r filename rest; do
        # Ignore empty or commented lines.
        if [[ -n $filename && ${filename:0:1} != "#" ]]; then
            echo "${filename}|${rest}" >> "${results_dir}/${filename}.emapper.annotations"
        fi
    done < $results_file
    # Add the info to each new annotation file and update the gene call id as required.
    for i in $results_dir/*.emapper.annotations; do
        if [[ $(basename $i) != "ALL_GENES.emapper.annotations" ]]; then
            sed -i "1e cat $tmp" $i # Insert the info at the start of the file
            sed -i 's/[^ 	]*|/g/' $i # Replace a string without spaces or tabs followed by |
        fi
    done
    rm -f "$tmp"
    echo "Done!"
    
    echo ""
    echo "  Beginning annotation import..."
    echo "--------------------------------------------------------------------------------"
    for i in ${in_dir}/*.db; do
        if [[ ! -e $i ]]; then
            echo "No contigs database files in the specified input directory!"
            exit 1
        fi
        filename=$(basename $i)
        echo "Processing $filename"
        anvi-script-run-eggnog-mapper --annotation "${results_dir}/${filename}.emapper.annotations" --use-version $emapper_version -c $i
    done
    echo "Done!"
fi

