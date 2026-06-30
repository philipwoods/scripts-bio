#!/usr/bin/env python

import os
import sys
import argparse
from Bio import SeqIO

masses = {
        'A': 89.09,
        'R': 174.20,
        'N': 132.12,
        'D': 133.10,
        'C': 121.15,
        'E': 147.13,
        'Q': 146.15,
        'G': 75.07,
        'H': 155.16,
        'I': 131.17,
        'L': 131.17,
        'K': 146.19,
        'M': 149.21,
        'F': 165.19,
        'P': 115.13,
        'O': 255.31,
        'U': 168.05,
        'S': 105.09,
        'T': 119.12,
        'W': 204.23,
        'Y': 181.19,
        'V': 117.15,
        'water': 18.02
        }

def main(args):
    output = []
    records = SeqIO.parse(args.fasta, 'fasta')
    for record in records:
        seq = str(record.seq).upper()
        bonds = len(seq) - 1
        mass = sum([masses[residue] for residue in seq]) - (masses['water'] * bonds)
        megadaltons = mass / 1000000
        format_str = "{id}\t{MDa:." + str(args.precision) + "f}"
        row = format_str.format(id=record.id, MDa=megadaltons)
        output.append(row)
    print('\n'.join(output), file=args.output)

if __name__ == "__main__":
    desc = ("Takes in a protein FASTA and reports the estimated mass of each protein in MDa.")
    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument('fasta', metavar='FASTA', help="FASTA file containing the sequence(s) of interest")
    parser.add_argument('-o', '--output', default=sys.stdout, help="Path to write output table. Default: STDOUT")
    parser.add_argument('--precision', default=2, type=int, help="Set the number of decimal places to report. Default: %(default)s")
    args = parser.parse_args()
    if not os.path.isfile(args.fasta):
        sys.exit("Specified file does not exist: {}".format(args.fasta))
    main(args)

