#!/usr/bin/env python

import os
import sys
import pandas as pd
import argparse
from collections import defaultdict

def parse_field(field):
    # Fields have the format X(N.NNN) where
    #   X is a single letter amino acid code
    #   N.NNN is a probability assigned to that amino acid
    parts = field.rstrip(")").split("(")
    return (parts[0], float(parts[1]))

def main(args):
    data = defaultdict(list)
    with open(args.rstnode, 'r') as f:
        # For each line, split into fields at whitespace
        for line in f:
            fields = line.split()
            # Record the position in the reconstructed sequence
            data['site'].append(fields[0])
            fields = fields[3:]
            # For the probability data, extract the residue and probability
            for field in fields:
                residue, prob = parse_field(field)
                data[residue].append(prob)
    # Convert the data dictionary into a dataframe, indexed on sequence position
    df = pd.DataFrame(data)
    df = df.set_index('site')
    # Add extra columns indicating the reconstructed sequence and the confidence in each residue
    df['rst'] = df.idxmax(axis=1)
    df['prob'] = df.max(axis=1, numeric_only=True)
    # Handle terminal or file output formatting
    if args.out is None:
        print(df.to_string(float_format='{:.3f}'.format))
    else:
        df.to_csv(args.out, sep='\t', float_format='{:.3f}'.format)

if __name__ == "__main__":
    desc = ("This script will convert the PAML rst output for a single node into a CSV "
            "format for ease of use in other programs. The input file should contain "
            "the PAML output table without any header lines.")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('rstnode', help="A text file containing the PAML rst output from one node")
    parser.add_argument('--out', '-o', help="A file to store CSV output in.")
    args = parser.parse_args()
    if not os.path.isfile(args.rstnode):
        sys.exit("File does not exist: {}".format(args.rstnode))
    main(args)

