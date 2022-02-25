#!/bin/bash

print_usage() {
    echo "THIS SCRIPT WORKS ON THE OUTPUT OF MY SCRIPT ANVI-RUN-KOFAMS.SH"
    echo "IF YOU USED THE BUILT-IN KEGG ANNOTATION FUNCTION IN ANVIO 7+ THIS WILL NOT WORK"
    echo ""
    echo "bash run-kegg-decoder.sh [-e EXT] -o OUT_DIR -i IN_DIR"
    echo "Looks in <in-dir> for KEGG output files created by anvi-run-kofams.sh"
    echo "The -e option lets you specify the file extension of the Kofam output. Default is ko"
}

print_version() {
    echo "Last modified 27 Jun 2020"
}

in_dir=""
out_dir=""
ext="ko"

while getopts ":hvi:o:e:" opt; do
    case $opt in
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        i)  in_dir="${OPTARG}"
            ;;
        o)  out_dir="${OPTARG}"
            ;;
        e)  ext="${OPTARG}"
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

echo ""
echo "  Beginning input file construction..."
echo "--------------------------------------------------------------------------------"

all_kos="${out_dir}/ALL_GENES.ko"

for i in ${in_dir}/*.${ext}; do
    if [[ ! -e $i ]]; then
        echo "No files with the given extension in the specified directory: $ext"
        exit 1
    fi
    filename=$(basename $i)
    # Remove extension from the file name.
    genomename=${filename%.*}
    # Replace all underscores with dashes.
    name=${genomename//_/-}
    # Put the edited name at the front of all lines in the ko file and output to the complete list.
    sed "s/^genecall_/${name}_/" $i >> $all_kos
    echo "Done adding $filename"
done

echo ""
echo "  Beginning KEGG decoding..."
echo "--------------------------------------------------------------------------------"

module load anaconda3
source activate keggdecoder

decoder_out="${out_dir}/KEGG-decoder.list"
KEGG-decoder -i $all_kos -o $decoder_out -v static
# Change the first column of the output to match the genome names rather than the ko file names.
sed -i 's/^kofams-//' $decoder_out
tmp=$(mktemp)
cp "$decoder_out" "$tmp" && gawk -v FS='\t' -v OFS='\t' '{gsub(/[-\.]/,"_",$1) ; print}' $tmp > $decoder_out
rm -f "$tmp"
# Change the column separators to conform to a tsv format.
sed -i 's/\t /\t/g' $decoder_out

echo "Done!"

