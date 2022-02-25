#!/bin/bash

print_usage() {
    echo "bash align-individual-marker-genes.sh [-c] [-t TREE_DIR -n NUM_THREADS -s SEED] -f MARKER_GENES -o OUT_DIR"
    echo "  The marker genes file should be generated from anvio without the --concatenate flag."
    echo "  The -c flag indicates that the input marker file is pangenome SCGs from"
    echo "  anvi-get-sequences-for-gene-clusters instead of markers from anvi-get-sequences-for-hmm-hits."
    echo "  By default, this program will only separate and align marker genes. If the -t option is used, trees will be"
    echo "  built for each marker using RAxML. The -n and -s options can be used to specify parameters to pass to RAxML."
    echo "  By default, 4 threads will be used, RAxML output will be placed in OUT_DIR, and a seed will be randomly generated."
}

print_version() {
    echo "Last modified 20 Jul 2020"
}

origin="hmms"
markers=""
out_dir="."
tree_dir=""
threads=4
seed=$RANDOM

while getopts ":hvct:n:s:f:o:" opt; do
    case $opt in
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        c)  origin="clusters"
            ;;
        t)  tree_dir="${OPTARG}"
            ;;
        n)  threads="${OPTARG}"
            ;;
        s)  seed="${OPTARG}"
            ;;
        f)  markers="${OPTARG}"
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

if [[ ! $CONDA_DEFAULT_ENV == "phylogeny" ]]; then
    echo "You must activate the phylogeny conda environment before running this."
    exit 1
fi

echo ""
echo "  Separating each marker gene into individual files..."
echo "--------------------------------------------------------------------------------"
python /data1/projects/pwoods/scripts/separate-marker-genes.py $markers $origin $out_dir

echo ""
echo "  Aligning each individual marker gene..."
echo "--------------------------------------------------------------------------------"
# Sequences taken from gene clusters are already aligned by anvio using muscle.
num_files=$(ls -1 ${out_dir}/marker-* | wc -l)
counter=0
if [[ $origin == "hmms" ]]; then
    echo "Looking in $out_dir for marker files."
    for f in $out_dir/marker-*; do
        tmp=$(mktemp)
        cp $f $tmp
        muscle -in $tmp -out $f -quiet
        rm -f $tmp
        ((counter++))
        base=$(basename $f)
        progress="($counter of $num_files)"
        printf "Completed aligning "
        printf "%-50s %10s\n" "$base" "$progress"
    done
fi
if [[ $origin == "clusters" ]]; then
    echo "Gene clusters exported from Anvio are already aligned."
fi

# If a tree output directory is provided, build trees.
if [[ -n $tree_dir ]]; then
    echo ""
    echo "  Building individual gene trees..."
    echo "--------------------------------------------------------------------------------"
    if [[ ! -d ${tree_dir} ]]; then
        mkdir ${tree_dir}
    fi
    counter=0
    real_tree_dir=$(realpath $tree_dir)
    echo $real_tree_dir
    for f in ${out_dir}/marker-*; do
        filename=$(basename $f)
        raxmlHPC-PTHREADS-AVX2 -T $threads -f a -p $seed -x $seed -# autoMRE -m PROTGAMMALG4X -s $f -n ${filename%.*} -w $real_tree_dir
        ((counter++))
        base=$(basename $f)
        progress="($counter of $num_files)"
        printf "Completed tree for "
        printf "%-50s %10s\n" "$base" "$progress"
    done
fi

echo "Done!"

