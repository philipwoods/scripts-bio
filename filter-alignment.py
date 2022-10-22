#!/usr/bin/env python

import os
import sys
import argparse
from Bio import AlignIO

def main(args):
    # Read in alignment file
    alignment = AlignIO.read(args.fasta, "fasta")
    ## Late argument validation
    # If no outgroup is given, the outgroup should be empty.
    if args.outgroup is None:
        args.outgroup = []
    out_set = set(args.outgroup)
    # If no ingroup is given, the ingroup should be everything not in the outgroup.
    if args.ingroup is None:
        in_set = set(range(len(alignment)))
        in_set = in_set.difference(out_set)
        args.ingroup = sorted(in_set)
    in_set = set(args.ingroup)
    if not in_set.isdisjoint(out_set):
        sys.exit("Sequences cannot be in both the ingroup and the outgroup.")
    # Handle list display argument
    if args.list:
        for i, record in enumerate(alignment):
            print("{index}\t{header}".format(index=i, header=record.id))
    print("Done")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', metavar='FASTA', help="FASTA file containing the alignment of interest")
    parser.add_argument('--list', '-l', action='store_true', help="Display indices of sequence headers in the FASTA")
    parser.add_argument('--ingroup', type=int, nargs='+', help="Specify indices of sequences to check for site conservation. Default: all")
    parser.add_argument('--outgroup', type=int, nargs='+', help="Specify indices of sequences which should not be conserved with the ingroup. Default: none")
    parser.add_argument('--identity', '-i', type=float, default=1, help="Set the minimum identity threshold for conservation. Default: %(default)s")
    args = parser.parse_args()
    # Argument validation
    if not os.path.isfile(args.fasta):
        sys.exit("Specified file does not exist: {}".format(args.fasta))
    if (args.identity < 0) or (args.identity > 1):
        sys.exit("Identity threshold must be a proportion between 0 and 1, inclusive.")
    main(args)

