#!/usr/bin/env python3
import sys
from collections import defaultdict
from Bio import SeqIO

def main(args):
    if "-h" in args or "--help" in args:
        print("Identifies identical sequence IDs in a FASTA file.")
        print("Takes one argument: input FASTA file.")
        sys.exit(1)

    duplication_record = defaultdict(int)
    for record in SeqIO.parse(args[0], "fasta"):
        duplication_record[str(record.id)] += 1

    print("Duplicate sequence IDs in " + args[0])
    print("Count\tSequence ID")
    for seqID, count in duplication_record.items():
        if count > 1:
            print("{0}\t{1}".format(count, seqID))

if __name__ == "__main__":
    main(sys.argv[1:])
