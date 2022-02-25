#!/usr/bin/env python3
import sys
from Bio import SeqIO

def main(args):
    if "-h" in args or "--help" in args:
        print("Truncates sequences in a FASTA file.")
        print("Takes three arguments: input FASTA file, start position, length.")
        sys.exit(1)

    start = int(args[1])
    end = int(args[2]) + start
    for record in SeqIO.parse(args[0], "fasta"):
        sliced = record[start:end]
        print(sliced.format("fasta"))

if __name__ == "__main__":
    main(sys.argv[1:])
