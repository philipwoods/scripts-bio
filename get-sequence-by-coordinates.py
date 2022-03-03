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
    # The provided start and stop positions are 1-indexed, but python is 0-indexed. Therefore we need to subtract 1 from each provided index.
    # The provided start and stop positions are inclusive, but python excludes the second slicing index. We need to add 1 to the second index.
    feature = contig[args.start-1:args.stop]
    if args.strand == '-':
        feature = feature.reverse_complement()
    # Handle sequence type and stop codon
    if args.nucleotide:
        if args.include_stop_codon is not None and not args.include_stop_codon:
            feature = feature[:-3]
    else:
        feature = feature.translate()
        if args.include_stop_codon is None or not args.include_stop_codon:
            feature = feature[:-1]
    # Set record information
    if len(args.name.split()) > 1:
        sys.exit("Provided name cannot contain whitespace: {}".format(args.name))
    feature.id = args.name
    feature.description = args.description
    # Write feature out
    out_path = os.path.join(args.output, args.name + ".fa")
    SeqIO.write(feature, out_path, 'fasta')

if __name__== "__main__":
    desc = ("This script will extract nucleotide or amino acid sequences from an input FASTA file based on "
            "a contig identifier, a strand direction, and start and stop positions within that contig. The "
            "resulting sequence will be written to a FASTA file with the provided feature name. The name "
            "and description will be used for the FASTA defline.")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('fasta', metavar='FASTA', help="A FASTA file containing contig sequences.")
    parser.add_argument('contig', help="The sequence identifier of the contig of interest.")
    parser.add_argument('start', type=int, help="The start position (1-indexed, inclusive) of the feature in the contig.")
    parser.add_argument('stop', type=int, help="The end position (1-indexed, inclusive) of the feature in the contig.")
    parser.add_argument('strand', metavar='strand', choices=['+','-'], help="The strand which the feature is on. Choices: {%(choices)s}")
    parser.add_argument('name', help="The name of the feature. Must not contain whitespace.")
    parser.add_argument('-o', '--output', default='.', help="The directory to store output in. Default: ./")
    parser.add_argument('-d', '--description', default='', help="Additional information about the feature e.g. functional or domain annotation(s). Remember to quote this if it contains whitespace.")
    parser.add_argument('-n', '--nucleotide', action='store_true', help="Output nucleotide sequence instead of amino acid sequence. Default: off")
    parser.add_argument('--include-stop-codon', type=int, choices=[1, 0], help="Choose whether to report the final stop codon in the output. Default: off (0) for AA sequences, on (1) for nt sequences.")
    args = parser.parse_args()
    if not os.path.isfile(args.fasta):
        sys.exit("File does not exist: {}".format(args.fasta))
    main(args)

