#!/usr/bin/env python3
import os
import sys
import argparse
from Bio import SeqIO

def main(args):
    for record in SeqIO.parse(args, 'fasta'):
        print(record.id, len(record), sep='\t')

if __name__ == "__main__":
    args = sys.argv[1]
    main(args)
