#!/bin/bash

print_usage() {
    echo "This function works properly versions of Anvio before v7. As of v7 and onwards"
    echo "it will not provide useful output"
    echo ""
    echo "Synopsis:"
    echo "  anvi-script-get-annotated-genome [-h] [-v] -c CONTIGS_DB -n GENOME_NAME"
    echo "                                   -f ANNOTATION_SOURCE -o OUT_DIR [-p PAN_DB]"
    echo ""
    echo "Description:"
    echo "  This script exports the contigs from the specified Anvi'o contigs database"
    echo "  as a fasta file, then constructs a feature annotation file in gff3 format"
    echo "  that corresponds to that fasta file. The annotation file includes the gene"
    echo "  caller ID and predicted product according to the given functional annotation"
    echo "  source. If an Anvi'o pangenome database is also provided, then the file will"
    echo "  also include the gene cluster ID. If a pangenome database is provided, then"
    echo "  the appropriate genome name must also be provided."
    echo ""
    echo "  The most useful annotation sources for genome gazing are probably Pfam,"
    echo "  COG20_FUNCTION, KOfam, EGGNOG_BACT, and Transfer_RNAs."
    echo "  For a more general overview, KEGG_Module, KEGG_Class, COG20_CATEGORY, and"
    echo "  COG20_PATHWAY may also be helpful."
    echo ""
    echo "Options:"
    echo "  -h: Display this help text and exit."
    echo "  -v: Display version text and exit."
    echo "  -p: Specify Anvi'o pangenome database to use."
    echo "  -n: Specify genome name to use."
    echo "  -c: Specify Anvi'o contigs database to use."
    echo "  -f: Specify annotation source for contigs database."
    echo "      To see supported sources, run anvi-export-functions -c CONTIGS_DB -l."
    echo "  -o: Specify output directory."
}

print_version() {
    echo "Last updated 22 October 2021"
}

pan_db=""
genome_name=""
contigs_db=""
annotation_src=""
out_dir=""

while getopts ":hvp:n:c:f:o:" opt; do
    case $opt in 
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        p)  pan_db="${OPTARG}"
            ;;
        n)  genome_name="${OPTARG}"
            ;;
        c)  contigs_db="${OPTARG}"
            ;;
        f)  annotation_src="${OPTARG}"
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

if [ -n "$pan_db" ] && [ ! -e "${pan_db}" ]; then
    echo "File does not exist: ${pan_db}"
    exit 1
fi
if [ -n "$pan_db" ] && [ -z "$genome_name" ]; then
    echo "You must provide a genome name to go with the pangenome."
    exit 1
elif [ -z "$pan_db" ] && [ -z "$genome_name" ]; then
    genome_name=$(basename $contigs_db .db)
fi
if [ ! -e "${contigs_db}" ]; then
    echo "File does not exist: ${contigs_db}"
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

genome_dir="${out_dir}/${genome_name}"
if [ ! -e "$genome_dir" ] || [ ! -d "$genome_dir" ]; then
    mkdir "$genome_dir"
fi
contigs_out="${genome_dir}/${genome_name}_contigs.fa"
functions_out="${genome_dir}/${genome_name}_functions_${annotation_src}.tab"
features_out="${genome_dir}/${genome_name}_features_${annotation_src}.gff3"

echo "Exporting contigs database..."
anvi-export-contigs -c "$contigs_db" -o "$contigs_out"
echo "Exporting feature annotations..."
anvi-get-sequences-for-gene-calls -c "$contigs_db" --export-gff3 -o "$features_out"
echo "Exporting functions..."
anvi-export-functions -c "$contigs_db" --annotation-sources "$annotation_src" -o "$functions_out"
if [ -n "$pan_db" ]; then
    pan_name=$(basename "$pan_db" .db)
    gc_out="${out_dir}/${pan_name}_gene_clusters.tab"
    filtered_gc_out="${genome_dir}/${genome_name}_gene_clusters.tab"
    echo "Exporting pangenome gene clusters..."
    anvi-export-table --table gene_clusters -o "$gc_out" "$pan_db"
    echo "Filtering gene clusters by genome name..."
    echo "Pangenome gene cluster file: $gc_out"
    echo "Search term: $genome_name"
    # Currently, the header for the gene clusters table is as follows
    # gene_caller_id	gene_cluster_id	genome_name	alignment_summary
    head -1 "$gc_out" | cut -f 1,2 > "$filtered_gc_out"
    grep -w "$genome_name" "$gc_out" | cut -f 1,2 >> "$filtered_gc_out"
    echo "Output file: $filtered_gc_out"
    echo "Adding information to feature annotation file..."
    python "${DIR}/add-info-to-gff3.py" "$features_out" "$functions_out" "$filtered_gc_out"
    rm -f "$gc_out"
else
    echo "Adding information to feature annotation file..."
    python "${DIR}/add-info-to-gff3.py" "$features_out" "$functions_out"
fi
exit
