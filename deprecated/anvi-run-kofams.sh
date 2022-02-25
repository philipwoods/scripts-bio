#!/bin/bash

print_usage() {
    echo "NOTE: THIS IS NOT NECESSARY WITH ANVIO 7 AND ABOVE."
    echo "USE THE BUILT-IN KEGG UTILITY WITH LATER VERSIONS OF ANVIO."
    echo "bash anvi-run-kofams.sh -c CONTIGS_DB -T NUM_THREADS -o OUT_DIR"
}

print_version() {
    echo "Last modified 9 Oct 2021"
}

contigs_db=""
out_dir=""
num_threads=""

while getopts ":hvc:o:t:" opt; do
    case $opt in
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        c)  contigs_db="${OPTARG}"
            ;;
        o)  out_dir="${OPTARG}"
            ;;
        t)  num_threads="${OPTARG}"
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

module load anaconda3
source activate anvio6

filename=$(basename ${contigs_db})

echo ""
echo "  Beginning gene call export..."
echo "--------------------------------------------------------------------------------"
# The parameter expansion removes the .db file extension.
gene_calls="${out_dir}/gene_calls_${filename%.*}.fa"
anvi-get-sequences-for-gene-calls -c $contigs_db -o $gene_calls --get-aa-sequences
sed -i 's/>/>genecall_/' $gene_calls

echo ""
echo "  Beginning KOfam annotation..."
echo "--------------------------------------------------------------------------------"
conda deactivate
module unload anaconda3
module load kofamscan

kofam_out="${out_dir}/kofams_${filename%.*}.ko"
exec_annotation --cpu=$num_threads -k /data1/sw/kofamscan/db/ko_list -p /data1/sw/kofamscan/db/profiles/ -f mapper -o $kofam_out $gene_calls

module unload kofamscan
module load anaconda3

anvi_table="${out_dir}/kofams_anvi_table_${filename%.*}.tsv"
python /data1/projects/pwoods/KEGG/GhostKoalaParser/KEGG-to-anvio --KeggDB /data1/projects/pwoods/KEGG/KO_Orthology_ko00001.txt -i $kofam_out -o $anvi_table
sed -i 's/KeggGhostKoala/KeggKofams/' $anvi_table

echo ""
echo "  Beginning annotation import..."
echo "--------------------------------------------------------------------------------"
source activate anvio6
anvi-import-functions -c $contigs_db -i $anvi_table
rm -f $gene_calls

echo ""
echo "Done!"

