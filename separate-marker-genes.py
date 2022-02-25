#!/usr/bin/env python3
import sys
import os.path
from collections import defaultdict
from Bio import SeqIO

def parse_defline(defline, form):
    """Takes in a defline string and extracts the unique marker name based on the format."""
    if form == "clusters":
        # Format is ENTRY#|gene_cluster:GC_########|genome_name:ORGANISM|gene_callers_id:####
        marker = defline.split("|")[1].split(":")[1]
        organism = defline.split("|")[2].split(":")[1]
        source = "SCG"
        other = defline.split("|")[3]
        return (marker, organism, source , other)
    if form == "hmms":
        # Format is MARKERNAME___HMMSOURCE___HASH bin_id:ORGANISM|source:HMMSOURCE|OTHERSTUFF
        marker = defline.split("___")[0]
        organism = defline.split("|")[0].split(":")[-1]
        source = defline.split("___")[1]
        other = "|".join(defline.split("|")[2:]) # Don't duplicate the source
        return (marker, organism, source, other)

def main(args):
    if "-h" in args or "--help" in args:
        print("Usage: python separate-marker-genes.py MARKERFILE ORIGIN OUT_DIR")
        print("This is a helper file for align-individual-sequences.sh.")
        print("ORIGIN should be 'hmms' or 'clusters' depending on the source of the marker file within anvio.")
        sys.exit(1)

    all_markers = args[0]
    origin = args[1]
    out_dir = args[2]

    if origin not in ["hmms", "clusters"]:
        print("Marker source not supported: %s" % origin)
        print("Use the --help flag for more information.")
        sys.exit(1)

    # Build a dictionary grouping all records in the input file by marker gene
    marker_records = defaultdict(list)
    for record in SeqIO.parse(all_markers, "fasta"):
        # The record ID is anything before a space in the defline, and the description is the whole defline.
        marker, organism, source, other = parse_defline(record.description, origin)
        # If the description doesn't include the ID, the defline will be "ID DESCRIPTION" in the exported file.
        record.id = organism
        if origin == "hmms":
            record.description = "source:%s|gene:%s|%s" % (source, marker, other)
        if origin == "clusters":
            record.description = "source:%s|gene_cluster:%s|%s" % (source, marker, other)
        key = "-".join(("marker", source, marker)) + ".fa"
        marker_records[key].append(record)
    # Go through each unique marker and export the corresponding records to a file.
    for filename in marker_records.keys():
        SeqIO.write(marker_records[filename], os.path.join(os.path.abspath(out_dir), filename), "fasta")

if __name__ == "__main__":
    main(sys.argv[1:])

