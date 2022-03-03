#!/usr/bin/env python3
import os
import sys
import argparse

def main(args):
    print(asdf)

if __name__== "__main__":
    desc = ("")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('fasta', metavar='FASTA', help="A FASTA file containing contig sequences.")
    parser.add_argument('contig', help="The sequence identifier of the contig of interest.")
    parser.add_argument('start', type=int, help="The start position of the feature in the contig.")
    parser.add_argument('stop', type=int, help="The end position of the feature in the contig.")
    parser.add_argument('name', help="The name of the feature.")
    parser.add_argument('-o', '--output', default='.', help="The directory to store output in. Default: ./")
    parser.add_argument('-n', '--nucleotide', action='store_true', help="Output nucleotide sequence instead of amino acid sequence. Default: off")
    args = parser.parse_args()
    if not os.path.isfile(args.fasta):
        sys.exit("File does not exist: {}".format(args.fasta))
    main(args)

