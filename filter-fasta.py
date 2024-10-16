#!/usr/bin/env python3
import sys
import os.path
import argparse
from Bio import SeqIO

def within_length_bounds(record, args):
    return (len(record) > args.min) and (len(record) < args.max)

def main(args):
    if args.list:
        for record in SeqIO.parse(args.fasta, "fasta"):
            print(record.id)
        sys.exit()
    if args.id_file is not None:
        with open(args.id_file) as f:
            file_ids = [line.strip() for line in f.readlines()]
        if args.id is None:
            args.id = file_ids
        else:
            args.id.extend(file_ids)
    output_records = []
    for record in SeqIO.parse(args.fasta, "fasta"):
        # The record ID is anything before a space in the defline
        if ((args.id is None) or (record.id in args.id)) and within_length_bounds(record, args):
            output_records.append(record)
    SeqIO.write(output_records, args.out, "fasta")

if __name__ == "__main__":
    desc = ("Takes in a multi-sequence FASTA file and searches it for the provided sequence IDs "
            "and/or by the provided minimum and maximum sequence lengths.")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("fasta", metavar="FASTA", help="A multi-sequence FASTA file.")
    parser.add_argument("--out", default=sys.stdout, help="A path to output the results. Default: stdout")
    parser.add_argument("--id", nargs='+', help="A list of sequence IDs to output. Default: output everything")
    parser.add_argument("--id-file", help="A file containing a list of sequence IDs, one per line, to output.")
    parser.add_argument("--list", "-l", action='store_true', help="List the available sequence IDs from the input.")
    parser.add_argument("--min", type=int, default=0, help="Minimum acceptable sequence length to output.")
    parser.add_argument("--max", type=int, default=100000000000000000, help="Maximum acceptable sequence length to output.")
    args = parser.parse_args()
    if not os.path.isfile(args.fasta):
        sys.exit("The specified file does not exist: {}".format(args.fasta))
    if args.id_file is not None and not os.path.isfile(args.id_file):
        sys.exit("The specified file does not exist: {}".format(args.id_file))
    main(args)

