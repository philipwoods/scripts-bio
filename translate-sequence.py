#!/usr/bin/env python3

import os
import sys
import argparse
from Bio import SeqIO

def is_nucleotide(sequence):
    # Define allowed characters for nucleotide sequences
    allowed = set(['A', 'a', 'T', 't', 'C', 'c', 'G', 'g', 'U', 'u'])
    # Create a set of all unique characters in the sequence
    chars = set(sequence)
    # Check whether all of the characters in the sequence are allowed nucleotides
    return chars.issubset(allowed)
    

def main(args):
    translations = []
    for record in SeqIO.parse(args.fasta, 'fasta'):
        if is_nucleotide(record.seq):
            translated = record.translate(
                id="translated_{}".format(record.id),
                description="translated_{}".format(record.description),
                to_stop=args.no_stop
            )
            translations.append(translated)
        else:
            translations.append(record)
    SeqIO.write(translations, args.output, 'fasta')

if __name__ == "__main__":
    desc = ""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('fasta', metavar='FASTA', help="A FASTA file containing nucleotide sequences.")
    parser.add_argument('--no_stop', action='store_false', help="Do not halt translation at stop codons but indicate them with '*'. Default: off")
    parser.add_argument('--output', '-o', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help="A file to store the translated output. Prints to stdout by default.")
    args = parser.parse_args()
    main(args)

