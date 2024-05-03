#!/usr/bin/env python3
import sys
import os.path
import argparse
from Bio import SeqIO

def main(args):
    records = SeqIO.parse(args.fasta, "fasta")
    # The record ID is anything before a space in the defline
    ids = [record.id for record in records]
    ids.sort()
    for seq_id in ids:
        # SeqIO.parse returns an iterator not a list, so we need to regenerate it
        for record in SeqIO.parse(args.fasta, "fasta"):
            if record.id == seq_id:
                SeqIO.write(record, sys.stdout, "fasta")

if __name__ == "__main__":
    desc = "Takes in a multi-sequence FASTA file and sorts the sequences by sequence ID."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("fasta", metavar="FASTA", help="A multi-sequence FASTA file.")
    args = parser.parse_args()
    if not os.path.isfile(args.fasta):
        sys.exit("The specified file does not exist: {}".format(args.fasta))
    main(args)

