#!/bin/bash

if [[ "$#" -ne 1 || $1 == "-h" || $1 == "--help" ]]; then
    echo "Usage: clean-gff-products [-h] <GFF>"
    echo ""
    echo "  Takes in GFF3 files produced by Anvi'o and escapes"
    echo "  reserved characters in the attribute fields. Edits"
    echo "  the input file in place."
    exit 1
fi

if [[ ! -f $1 ]]; then
    echo "File does not exist: $1"
    exit 1
fi

# Escape semicolons that are followed by anything other than certain field names.
# Need to use perl because sed and awk don't have lookahead.
perl -pi -e 's/;(?!ID|Name|Alias|db_xref|[Pp]roduct|Note)/ %3B /g' $1
# Escape all commas
sed -i 's/,/ %2C /g' $1

exit
