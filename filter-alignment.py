#!/usr/bin/env python

import os
import sys
import argparse
from Bio import AlignIO

def main(args):
    alignment = AlignIO.read(args.fasta, "fasta")
    if args.list:
        for i, record in enumerate(alignment):
            print("{index}\t{header}".format(index=i, header=record.id))
    print("Done")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', metavar='FASTA', help="FASTA file containing the alignment of interest")
    parser.add_argument('--list', '-l', action='store_true', help="Display indices of sequence headers in the FASTA")
    parser.add_argument('--identity', '-i', type=float, default=1, help="Set the minimum identity threshold for conservation. Default: %(default)s")
    args = parser.parse_args()
    if not os.path.isfile(args.fasta):
        sys.exit("Specified file does not exist: {}".format(args.fasta))
    if (args.identity < 0) or (args.identity > 1):
        sys.exit("Identity threshold must be a proportion between 0 and 1, inclusive.")
    main(args)

