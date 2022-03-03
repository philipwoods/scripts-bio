#!/usr/bin/env python3
import os
import sys
import argparse
from Bio import SeqIO

def main(args):
    # Identify the contig of interest
    contig = None
    for record in SeqIO.parse(args.fasta, 'fasta'):
        if record.id == args.contig:
            contig = record
            break
    if contig is None: # If we didn't find the contig, we can't proceed
        sys.exit("The specified contig ID was not found: {}".format(args.contig))
    # Get the feature out of the contig
    # The first index is inclusive, but the second is not
    # We want both to be inclusive
    feature = contig[start:stop+1]
    if args.strand == '-':
        feature = feature.reverse_complement()
    if not args.nucleotide:
        feature = feature.translate()
    # Set record information
    if len(args.name.split()) > 1:
        sys.exit("Provided name cannot contain whitespace: {}".format(args.name))
    feature.id = args.name
    feature.description = args.description
    # Write feature out
    out_path = os.path.join(args.output, args.name + ".fa")
    SeqIO.write(feature, out_path, 'fasta')

if __name__== "__main__":
    desc = ("")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('fasta', metavar='FASTA', help="A FASTA file containing contig sequences.")
    parser.add_argument('contig', help="The sequence identifier of the contig of interest.")
    parser.add_argument('start', type=int, help="The start position of the feature in the contig.")
    parser.add_argument('stop', type=int, help="The end position of the feature in the contig.")
    parser.add_argument('strand', metavar='strand', choices=['+','-'], help="The strand which the feature is on. Choices: {%(choices)s}")
    parser.add_argument('name', help="The name of the feature. Must not contain whitespace.")
    parser.add_argument('-o', '--output', default='.', help="The directory to store output in. Default: ./")
    parser.add_argument('-d', '--description', default='', help="Additional information about the feature e.g. functional or domain annotation(s).")
    parser.add_argument('-n', '--nucleotide', action='store_true', help="Output nucleotide sequence instead of amino acid sequence. Default: off")
    args = parser.parse_args()
    if not os.path.isfile(args.fasta):
        sys.exit("File does not exist: {}".format(args.fasta))
    main(args)

