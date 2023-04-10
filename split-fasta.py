#!/usr/bin/env python3
import sys
import os.path
import argparse
from Bio import SeqIO

def main(args):
    if args.id_file is not None:
        with open(args.id_file) as f:
            file_ids = [line.strip() for line in f.readlines()]
        if args.id is None:
            args.id = file_ids
        else:
            args.id.extend(file_ids)
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
    parser.add_argument("--id-file", help="A file containing a list of sequence IDs, one per line, to output")
    args = parser.parse_args()
    if not os.path.isfile(args.fasta):
        sys.exit("The specified file does not exist: {}".format(args.fasta))
    if args.id_file is not None and not os.path.isfile(args.id_file):
        sys.exit("The specified file does not exist: {}".format(args.in_file))
    main(args)

