#!/bin/bash

usage() {
    echo "finish-alphafold.sh -d DIR"
    echo ""
    echo "Takes in the output directory created by alphafold2 which"
    echo "contains the model files that it generated. Renames these"
    echo "so they are easier to work with in PyMOL."
}

dir=""

while getopts ":hd:" opt; do
    case $opt in 
        h)  usage
            exit 1
            ;;
        d)  dir=${OPTARG}
            ;;
        \?) echo "Invalid option: -${OPTARG}. Use the -h option for more information.">&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument. Use the -h option for more information.">&2
            exit 1
            ;;
    esac
done

if [[ ! -d "$dir" ]]; then
    echo "Directory does not exist: ${dir}">&2
    exit 1
fi

# We will assume the directory name is a useful description of the contents e.g. the organism it's from
description=$(basename $dir)
for f in $dir/*.pdb; do
    # In raw Alphafold output, this will be either ranked_# or relaxed_# or unrelaxed_#
    filename=$(basename $f)
    # If the filename has the description (directory name) in it, we don't want to add it again.
    if [[ ! $filename =~ $description ]]; then
        mv $f "${dir}/${description}_${filename}"
    fi
done

