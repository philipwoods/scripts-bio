#!/bin/bash

print_usage() {
    echo ""
    echo "Synopsis:"
    echo "  find-hemec-beta-barrels.sh [-h] [-v] [-s] -c CONTIGSDB -o OUTDIR -a ANNOTATION-SOURCES [-r RELATIVES-DIR]"
    echo ""
    echo "Description:"
    echo "  This script identifies potential heme binding proteins in the input CONTIGSDB. It"
    echo "  summarizes the distribution of heme c binding sites, provides the protein sequences"
    echo "  of potential heme binding proteins, and provides sequences for a few proteins adjacent"
    echo "  to each of these potential heme proteins. It also uses the provided ANNOTATION-SOURCES"
    echo "  to provide information about the potential function of the identified proteins."
    echo "  If a directory is provided containing genomic sequences or contigs databases for related"
    echo "  organisms, summary information about the distribution of heme binding motifs in those"
    echo "  organisms will also be computed for comparison."
    echo ""
    echo "  By default, this search uses a loose definition for the heme binding motif. To restrict"
    echo "  the search to only the canonical CXXCH, use the strict flag -s."
    echo ""
    echo "  Make sure to only use annotation sources with exactly one record per gene such as"
    echo "  COG_FUNCTION, EGGNOG_BACT, or KeggKofams. Pfam is not a suitable annotation source."
    echo "  To check what sources are available, run anvi-export-functions -l -c CONTIGSDB."
    echo "  To use multiple sources, provide a comma separated list i.e. source1,source2,etc."
    echo ""
    echo "Options:"
    echo "  -h: Display this usage information and exit."
    echo "  -v: Display version information and exit."
    echo "  -s: Search for hemes using a strict definition of the heme binding motif."
    echo "  -c: Specify Anvi'o contigs database for the organism of interest."
    echo "  -o: Specify an existing directory to place output files in."
    echo "  -a: Specify one or more annotation sources to use for functional information."
    echo "  -r: Specify a directory containing genomic fasta files or contigs databases for related"
    echo "      organisms. These can be within nested directories e.g. as in NCBI data archives."
    echo ""
}

print_version() {
    echo "Last updated 15 September 2020"
}

summarize_counts() {
    local count=$1
    local summary=$2
    local histogram=$3
    sed 1d "$count" | cut -f 2 | sort -g | uniq -c | sed 's/^ *//' | sed 's/\([0-9]*\) \([0-9]*\)/\2\t\1/' > "$summary"
    echo -e "Heme count\tFrequency" > "$histogram"
    perl -lane 'print $F[0], "\t", "=" x $F[1]' "$summary" >> "$histogram"
    local total=$(awk '{s+=$2} END {print s}' $summary)
    echo -e "Total\t$total" >> $summary
    sed -i "1i Heme count\tFrequency" "$summary"
}

contigsdb=""
outdir=""
annotations=""
relatives=""
motif="C.{2,5}CH|C..CK"

while getopts ":hvc:o:a:r:s" opt; do
    case $opt in
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        c)  contigsdb="${OPTARG}"
            ;;
        o)  outdir=$(realpath ${OPTARG})
            ;;
        a)  annotations="${OPTARG}"
            ;;
        r)  relatives=$(realpath ${OPTARG})
            ;;
        s)  motif="C..CH"
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
    echo "You need to activate an anvio conda environment before running this."
    exit 1
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

genome=$(basename $contigsdb .db)

# Get protein sequences for all genes in the genome
echo "Exporting protein sequences and functions..."
allgenesfile="${outdir}/${genome}_ALLGENES.faa"
anvi-get-sequences-for-gene-calls -c $contigsdb --get-aa-sequences -o "$allgenesfile"
funcsbase="${outdir}/${genome}_functions"
for fn in $(echo $annotations | tr "," " "); do
    funcs="${funcsbase}_${fn}.txt"
    anvi-export-functions -c $contigsdb --annotation-sources $fn -o "$funcs"
    tmpfuncs=$(mktemp)
    cut -f 2 --complement $funcs | sed 's/\t/|/g' | sed '1d' > $tmpfuncs
    mv $tmpfuncs $funcs
done

# Find the sequences of the genes with these motifs
echo "Identifying protein sequences with heme c binding motifs..."
searchfile="${outdir}/${genome}_heme_proteins.faa"
echo -n "" > "$searchfile"
$DIR/prot_motif_search.pl "$allgenesfile" "$motif" "$searchfile"
for fn in $(echo $annotations | tr "," " "); do
    searchfunc="${searchfile%.*}_${fn}.faa"
    funcs="${funcsbase}_${fn}.txt"
    cp "$searchfile" "${searchfunc}"
    python $DIR/modify_fasta_defline.py "$searchfunc" "$genome" "$funcs"
    rm -f $funcs
done

# Find the frequency of heme binding motifs in these genes and produce summary files
echo "Finding distribution of heme c binding motifs in genome of interest..."
countsfile="${outdir}/${genome}_counts.tsv"
summaryfile="${outdir}/${genome}_counts_summary.tsv"
histogramfile="${outdir}/${genome}_counts_histogram.tsv"
echo -n "" > "$countsfile"
$DIR/prot_motif_count.pl "$allgenesfile" "$motif" "$countsfile"
summarize_counts "$countsfile" "$summaryfile" "$histogramfile"

# If relatives are provided, do the same summary for them
if [[ -n $relatives ]]; then
    echo "Finding distribution of binding motifs in relatives..."
    relativesummary="${outdir}/${genome}_relatives_summary.tsv"
    echo -n "" > "$relativesummary"
    for i in $(find $relatives -name "*.fna" -o -name "*.db"); do
        base=${i%.*} # This includes the whole path except the file extension
        ext=${i##*.} # Get the extension i.e. fna or db, not .fna or .db
        reldb="${base}.db"
        if [[ $ext == "fna" ]] && [[ -e $reldb ]]; then
            continue
        elif [[ $ext == "fna" ]] && [[ ! -e $reldb ]]; then
            reformatted="${i}.reformatted"
            anvi-script-reformat-fasta -o $reformatted --simplify-names $i
            anvi-gen-contigs-database -f "$reformatted" -n "heme_comparison" -o "$reldb"
        fi
        seqs="${base}_ALLGENES.faa"
        count="${base}_counts.tsv"
        summary="${base}_counts_summary.tsv"
        hist="${base}_counts_histogram.tsv"
        echo -n "" > "$count"
        anvi-get-sequences-for-gene-calls -c "$reldb" --get-aa-sequences -o "$seqs"
        $DIR/prot_motif_count.pl "$seqs" "$motif" "$count"
        summarize_counts "$count" "$summary" "$hist"
        cat $summary >> "$relativesummary"
        echo "" >> "$relativesummary"
        rm -f $seqs
    done
fi

# Get info about the local surroundings of these genes
echo "Exporting gene loci..."
locibase="${outdir}/${genome}_heme_loci"
for fn in $(echo $annotations | tr "," " "); do
    echo -n "" > "${locibase}_${fn}.faa"
done
geneids=$(grep \> "$countsfile" | sed 's/>//' | cut -f 1 | tr '\n' ',' | sed 's/,$//')
locusdir="${outdir}/heme_loci/"
if [[ ! -e $locusdir ]] || [[ ! -d $locusdir ]]; then
    mkdir "$locusdir"
fi
anvi-export-locus -c $contigsdb --gene-caller-ids "$geneids" -o "$locusdir" -O heme_locus -n 2,2 -W
echo "Collecting protein information from each locus..."
for i in ${locusdir}/*.db; do
    locus=$(basename $i .db)
    proteinfile="${locusdir}/${locus}_proteins.faa"
    anvi-get-sequences-for-gene-calls --get-aa-sequences -c $i -o "$proteinfile"
    locusfuncbase="${locusdir}/${locus}_functions"
    for fn in $(echo $annotations | tr "," " "); do
        locusfuncs="${locusfuncbase}_${fn}.txt"
        anvi-export-functions -c $i --annotation-sources $fn -o "$locusfuncs"
        tmp=$(mktemp)
        cut -f 2 --complement $locusfuncs | sed 's/\t/|/g' | sed '1d' > $tmp
        mv $tmp $locusfuncs
        locusproteins="${proteinfile%.*}_${fn}.faa"
        cp $proteinfile $locusproteins
        python $DIR/modify_fasta_defline.py "$locusproteins" "${locus}" "$locusfuncs"
        cat "$locusproteins" >> "${locibase}_${fn}.faa"
        echo "" >> "${locibase}_${fn}.faa"
        rm -f "$locusfuncs"
    done
    rm -f $proteinfile
done

rm -f $searchfile
rm -f $allgenesfile

# Closing remarks
echo ""
echo "------------------------------------------------------------"
echo "  To identify potential transmembrane beta barrel domains,"
echo "  paste the sequences into the PRED-TMBB2 tool, found at:"
echo "    http://www.compgen.org/tools/PRED-TMBB2"
echo "------------------------------------------------------------"
echo ""

