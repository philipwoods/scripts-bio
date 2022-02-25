#!/usr/bin/env python3
import sys
import os.path
import argparse
from Bio import SeqIO

def main(args):
    for record in SeqIO.parse(args.fasta, "fasta"):
        # The record ID is anything before a space in the defline
        if (args.id is None) or (record.id in args.id):
            filename = "{}.fa".format(record.id)
            SeqIO.write(record, os.path.join(os.path.abspath(args.dir), filename), "fasta")

if __name__ == "__main__":
    desc = ("Takes in a multi-sequence FASTA file and splits it into one or more single-sequence "
            "FASTA files. The resulting files are named according to their sequence IDs. ")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("fasta", metavar="FASTA", help="A multi-sequence FASTA file.")
    parser.add_argument("-d", "--dir", default="", help="A directory to use for the output. Default: current directory")
    parser.add_argument("--id", nargs='+', help="A list of sequence IDs to output. Default: output everything")
    args = parser.parse_args()
    if not os.path.isfile(args.fasta):
        sys.exit("The specified file does not exist: {}".format(args.fasta))
    main(args)

