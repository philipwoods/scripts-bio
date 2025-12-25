#!/usr/bin/env python3
import sys
import os.path
import argparse
import pandas as pd
from Bio import SeqIO

def main(args):
    if args.list:
        for record in SeqIO.parse(args.fasta, "fasta"):
            print(record.id)
        sys.exit()
    name_df = pd.read_table(args.correspondence)
    if list(name_df.columns) != ['old_name', 'new_name']:
        sys.exit("The correspondence file must have two columns named 'old_name' and 'new_name'")
    output_records = []
    for record in SeqIO.parse(args.fasta, "fasta"):
        # The record ID is anything before a space in the defline
        if record.id not in name_df['old_name'].array:
            print("Sequence ID not found: {}".format(record.id))
            continue
        new_name = name_df[name_df['old_name']==record.id]['new_name'].to_numpy()[0]
        record.id = new_name
        output_records.append(record)
    SeqIO.write(output_records, args.out, "fasta")

if __name__ == "__main__":
    desc = ("Takes in a multi-sequence FASTA file and searches it for the provided sequence IDs "
            "and/or by the provided minimum and maximum sequence lengths.")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("fasta", metavar="FASTA", help="A multi-sequence FASTA file.")
    parser.add_argument("correspondence", help="A multi-sequence FASTA file.")
    parser.add_argument("--out", default=sys.stdout, help="A path to output the results. Default: stdout")
    parser.add_argument("--list", action='store_true', help="List sequence IDs in the provided fasta file.")
    args = parser.parse_args()
    if not os.path.isfile(args.fasta):
        sys.exit("The specified file does not exist: {}".format(args.fasta))
    if not os.path.isfile(args.correspondence):
        sys.exit("The specified file does not exist: {}".format(args.correspondence))
    main(args)

