#!/bin/bash

# TO DO maybe condense the process so this also includes the call to anvi-analyze-synteny?

print_usage() {
    echo "Synopsis:"
    echo "  anvi-script-extract-synteny [-h] [-v] [-d] -i LOCI -o OUT_DIR -n PREFIX"
    echo ""
    echo "Description:"
    echo "  This script takes the output of anvi-analyze-synteny and makes it informative."
    echo "  It searches that output for synteny groups that are present in multiple genomes"
    echo "  and creates a summary table of them, indicating the contents of the group, the"
    echo "  length of the group, and the genomes which contain it, as well as the number"
    echo "  of times it is present in each of those genomes."
    echo ""
    echo "  IMPORTANT: USE OF THE -d FLAG GREATLY INCREASES COMPUTATION TIME"
    echo ""
    echo "  The -d flag changes the output in an important way. If the -d flag is provided,"
    echo "  any subsets of a synteny record within a genome will be removed. In other words,"
    echo "  when this flag is provided, any record present in the output will represent a"
    echo "  fully unique set of genes within that genome. When this flag is off, many of the"
    echo "  records for a given genome may be subsets of another longer record. Deduplicating"
    echo "  the records may be useful by paring down the amount of data in the output, but"
    echo "  leaving all records including subsets may make it easier to notice when organisms"
    echo "  have regions of shared synteny that differ by a few genes at the end."
    echo ""
    echo "Options:"
    echo "  -h: Display this help text and exit."
    echo "  -v: Display version text and exit."
    echo "  -d: Deduplicate the synteny records in the output."
    echo "  -i: Specify the input (output from anvi-analyze-synteny)"
    echo "  -o: Specify output directory."
    echo "  -n: Specify a prefix for the summary output file."
}

print_version() {
    echo "Last updated 13 November 2021"
}

infile=""
outdir=""
prefix=""
dedup=""

while getopts ":hvdi:o:n:" opt; do
    case $opt in 
        h)  print_usage
            exit 1
            ;;
        v)  print_version
            exit 1
            ;;
        d)  dedup="-deduplicated"
            ;;
        i)  infile="${OPTARG}"
            ;;
        o)  outdir="${OPTARG}"
            ;;
        n)  prefix="${OPTARG}"
            ;;
        \?) echo "Invalid option: -${OPTARG}. Use the -h option for more information.">&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument. Use the -h option for more information.">&2
            exit 1
            ;;
    esac
done

if [[ ! -e "$infile" ]]; then
    echo "File does not exist: ${infile}"
    exit 1
fi
if [[ ! -d "$outdir" ]]; then
    echo "Directory does not exist: ${outdir}"
    exit 1
fi
if [[ -z "$prefix" ]]; then
    echo "Please provide an informative prefix for the output file."
    exit 1
fi

# The input file (as of anvio v7.1) has the following format:
# ngram    count    contig_db_name    N    number_of_loci
# Each line holds a record of a set of ordered genes in a single genome.
# The ngram is the description of the genes in order.
# The count is the number of times it occurs in the specified genome.
# The contig_db_name is the genome this ngram occurs in.
# N is the length of the ngram (how many genes in order).
# The number_of_loci is the number of genomes being compared.

echo "Please be patient! Each step can take quite some time."
outfile="${outdir}/${prefix}-synteny-condensed${dedup}.tsv"

## The first step is to filter out any ngrams that don't occur in multiple genomes.
## If the ngram is only present in one genome, we're not interested.
echo "Identifying common synteny across genomes..."
tmp1=$(mktemp)
tmp2=$(mktemp)
# The first command sorts so things with identical first columns (i.e. synteny) are adjacent. This is important for the use of uniq later.
# The second command moves the first field to the end of the line (and changes the field separator to a space).
# The third command removes the empty leading field and whitespace created by the second command.
# The fourth command picks out all lines which are duplicates beyond the fourth field and groups them,
#   leaving an empty line between groups. The last field (i.e. the synteny) is the fifth field.
sort "${infile}" | awk '{first = $1; $1=""; print $0, first}' | awk '{$1=$1};1' | uniq -f4 --all-repeated=separate > $tmp1
# The first command puts records into "synteny_id   length   genome_name   num_occurrences" format,
#   but also introduces empty fields separated by tabs to the blank lines.
# The second command clears the whitespace from formerly empty lines and sets the field separator to tab.
awk '{last = $NF; print last, $3, $2, $1}' $tmp1 | awk 'OFS="\t" {$1=$1};1' > $tmp2
rm -f $tmp1

## Once we have the 'singleton' ngrams removed, we can check for duplicates if necessary.
if [[ -n "${dedup}" ]]; then # If the dedup flag was provided...
    echo "Deduplicating synteny records..."
    # Get real directory where script file is located,
    # even if the script was called via a symlink
    SOURCE=${BASH_SOURCE[0]}
    while [ -L "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
      DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )
      SOURCE=$(readlink "$SOURCE")
      [[ $SOURCE != /* ]] && SOURCE=$DIR/$SOURCE # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
    done
    DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )
    tmp3=$(mktemp)
    # Pass tmp2 as the input and tmp3 as a location to write the output.
    python "${DIR}/deduplicate-synteny.py" ${tmp2} ${tmp3}
    # The output of this python script produces a file with many lines in a row with only whitespace.
    # We can use the following to remove the whitespace and extra lines.
    awk 'OFS="\t" {$1=$1};1' "${tmp3}" | cat -s > "${tmp2}"
fi

## After removing singleton ngrams (and maybe deduplicating), we can group the genomes with common synteny together.
echo "Separating synteny groups..."
# Prints each group to a new file (e.g. group-1.txt, group-2.txt,...)
awk -v dir="${outdir}" -v RS= '{print > (dir "/group-" NR ".txt")}' $tmp2
rm -f $tmp2

echo "Compiling output..."
# Print a header to the output file
echo "synteny	length	num_genomes	genome_names(num_occurrences)" > $outfile
# Print one line per synteny group to the output file
for groupfile in ${outdir}/group-*.txt; do
    # Extract the synteny ID for this group (same on all lines)
    synteny=$(cut -f 1 "${groupfile}" | uniq)
    # Extract the length of the group (same on all lines)
    length=$(cut -f 2 "${groupfile}" | uniq)
    # Count the number of organisms with this synteny group
    count=$(cut -f 1 "${groupfile}" | uniq -c | cut -f 1)
    # Print the relevant fields to the output file
    printf $synteny >> $outfile
    printf '\t' >> $outfile
    printf $length >> $outfile
    printf '\t' >> $outfile
    printf $count >> $outfile
    printf '\t' >> $outfile
    # Format the last column of the output properly.
    #   First, cut out the columns for genome names and the occurrences per genome.
    #   Format this information like "genome-name(num-occurrences)" on each line.
    #   Use tr to condense the multiple lines into a comma-separated list.
    #   Finally, replace the trailing comma with a newline.
    cut -f 3,4 "${groupfile}" | sed 's/\t/(/' | sed 's/$/)/' | tr "\n" "," | sed 's/,$/\n/' >> $outfile
    rm -f "${groupfile}"
done

echo "Done!"

exit
