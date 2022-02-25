#!/usr/bin/env python3
import sys
import os
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# For the record, the defline formats for anvio output from anvi-get-sequences-for-* are:
# hmms
# >GENOME_NAME num_genes:##|genes:G1,G2,G3,...|separator:\X\X\X
# gene-clusters
# >GENOME_NAME gene_clusters:GC1,GC2,GC3,...|separator:\X\X\X
# gene-calls
# >GENE_CALLER_ID contig:c_######;start:####;stop:####;direction:r;rev_compd:True;length:###
#
# The output of my script separate-marker-genes.py formats deflines like so:
# hmms
# >GENOME_NAME source:HMM_SOURCE|gene:G1|OTHER_STUFF
# clusters
# >GENOME_NAME source:HMM_SOURCE|gene_cluster:GC1|OTHER_STUFF
#
# The output of my script align-genes-by-id.py formats deflines like so:
# >GENOME_NAME gene:G1|gene_callers_id:##

def main(args):
    # Handle multiple possible input sources
    files = []
    if args.gene_directory:
        gene_dir = os.path.abspath(args.gene_directory)
        files = [os.path.join(gene_dir, f) for f in os.listdir(gene_dir) if (os.path.splitext(f)[-1] == args.file_extension)]
        if len(files) < 2:
            sys.exit("There were too few files with the extension {x} to concatenate in the specified directory:\n{d}".format(x=args.file_extension, d=args.gene_directory))
    if args.gene_files:
        files = [os.path.abspath(f) for f in args.gene_files]
    # Create a regex to find the start of the gene names list in the deflines.
    # Matches gene: genes: gene_cluster: gene_clusters: preceded by | or whitespace
    gene_list_pattern = re.compile(r"(?<=[\s\|])gene(s?|_cluster(s?)):")
    # Import all alignments from all organisms and get information about their lengths
    old_records = {}   # Will have the format {GENOME_NAME: {GENE_NAME: record}}
    alignment_lengths = {}  # Will have the format {GENE_NAME: length}
    for f in files:
        for record in SeqIO.parse(f, "fasta"):
            # Look for the gene list in the defline
            match = gene_list_pattern.search(record.description)
            # Exit if it's not there
            if match is None:
                sys.exit("Error: Missing gene name(s) for {genome}".format(genome=record.id))
            # If it is there, get the gene list and use it to store the record and its length
            index = match.start()
            alignment_name = record.description[index:].split("|")[0].split(":")[-1]
            # All alignments of a gene should be the same length, so we only need to add it once
            if alignment_name not in alignment_lengths:
                alignment_lengths[alignment_name] = len(record.seq)
            # Add the record to the appropriate dictionary
            if record.id in old_records:
                old_records[record.id][alignment_name] = record
            else: # Need to add the nested dictionary if it's not already present
                old_records[record.id] = {alignment_name: record}
    # Get the complete list of alignments
    alignments = alignment_lengths.keys()
    # Prepare to create a new record object for each organism
    new_records = []
    for genome_name in old_records:
        genome_sequences = []
        for alignment in alignments:
            if alignment in old_records[genome_name]:
                genome_sequences.append(old_records[genome_name][alignment].seq)
            else:
                genome_sequences.append("-" * alignment_lengths[alignment])
        # Concatenate gene alignments using the separator
        sep = Seq(args.separator)
        new_seq = sep.join(genome_sequences)
        # Construct new description and SeqRecord object
        gene_str = ",".join(alignments)
        new_desc = "num_genes:{count}|genes:{genes}|separator:{separator}".format(count=len(alignments), genes=gene_str,separator=args.separator)
        new_record = SeqRecord(new_seq, id=genome_name, description=new_desc)
        new_records.append(new_record)
    # Export the concatenated records
    SeqIO.write(new_records, os.path.abspath(args.output), "fasta")

if __name__ == "__main__":
    # Handle arguments
    func_desc = ("Takes in a directory containing protein sequence alignment files and "
                 "concatenates the alignments. The alignment deflines must use gene: or "
                 "genes: or gene_cluster: or gene_clusters: to specify the gene name(s).")
    parser = argparse.ArgumentParser(description=func_desc, allow_abbrev=False)
    parser.add_argument('-s', '--separator', default="XXX", help="The sequence to use to separate the alignments in the concatenation. Default is XXX")
    parser.add_argument('-x', '--file-extension', default=".fa", help="The file extension of the alignments. Only used with the -g argument. Other files will be ignored. Default is .fa")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-d', '--gene-directory', help="A directory containing protein alignment files.")
    group.add_argument('-f', '--gene-files', nargs='+', help="One or more protein alignment files.")
    parser.add_argument('-o', '--output', required=True, help="Specify the concatenated alignment output file.")
    args = parser.parse_args()
    if args.gene_directory and not os.path.isdir(args.gene_directory):
        sys.exit("The specified directory does not exist: {}".format(os.path.abspath(args.gene_directory)))
    if args.gene_files and (len(args.gene_files) < 2):
        sys.exit("You must provide multiple alignment files to concatenate.")

    main(args)

