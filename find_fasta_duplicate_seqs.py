#!/usr/bin/env python3
import sys
from collections import defaultdict
from Bio import SeqIO

def main(args):
    if "-h" in args or "--help" in args:
        print("Identifies identical sequences in a FASTA file.")
        print("Takes one argument: input FASTA file.")
        sys.exit(1)

    duplication_record = defaultdict(list)
    for record in SeqIO.parse(args[0], "fasta"):
        duplication_record[str(record.seq)].append(record.id)

    print("Duplicate sequences in " + args[0])
    print("Grouped lines are FASTA record IDs of identical sequences.")
    for entry in duplication_record.values():
        if len(entry) > 1:
            print("")
            print("\n".join(entry))

if __name__ == "__main__":
    main(sys.argv[1:])
