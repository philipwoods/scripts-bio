#!/usr/bin/env python3
import sys
import os.path
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(args):
    output_records = []
    for record in SeqIO.parse(args.fasta, "fasta"):
        # The record ID is anything before a space in the defline
        new_id = record.id + "_revcomp"
        output_records.append(record.reverse_complement(id=new_id, description=False))
    SeqIO.write(output_records, args.output, "fasta")

if __name__ == "__main__":
    desc = ("Takes in a multi-sequence FASTA file and outputs the reverse complement "
            "of each sequence, appending '_revcomp' to the sequence ID.")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("fasta", metavar="FASTA", help="A multi-sequence FASTA file.")
    parser.add_argument("output", help="Path to use when writing output file.")
    args = parser.parse_args()
    if not os.path.isfile(args.fasta):
        sys.exit("The specified file does not exist: {}".format(args.fasta))
    main(args)

