#!/bin/bash

usage() {
    echo "finish-tmhmm.sh [-q DPI] -d DIR"
    echo ""
    echo "Takes in the output directory created by tmhmm which contains"
    echo "the .plp and .gnuplot files that it generated. Uses these to"
    echo "create simple figures at the provided DPI (default 72)."
}

dir=""
quality=72

while getopts ":hd:q:" opt; do
    case $opt in 
        h)  usage
            exit 1
            ;;
        d)  dir=${OPTARG}
            ;;
        q)  quality=${OPTARG}
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

for f in $dir/*.gnuplot; do
    name=${f%.gnuplot}
    gnuplot $f # Creates an encapsulated postscript file ${f%.gnuplot}.eps
    convert -density $quality ${name}.eps ${name}.gif
done

