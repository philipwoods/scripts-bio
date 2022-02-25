#!/usr/bin/env python3
import os
import sys
import argparse

def main(args):
    output = []
    p_col = 2
    alpha = float(args.alpha)
    with open(args.input, 'r') as f:
        if not args.noheader:
            output.append(f.readline())
        for line in f:
            linelist = line.rstrip().split('\t')
            if float(linelist[p_col]) < alpha:
                output.append(line)
    print("".join(output))

if __name__ == "__main__":
    # Handle arguments
    desc = "Filters a coinfinder output file by a p-value threshold."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("alpha", help="The maximum p-value to include in the output.")
    parser.add_argument("input", help="The *_pairs.tsv output file from coinfinder.")
    parser.add_argument("--noheader", action='store_true', help="Specify that the input has no header row.")
    # Maybe change these to accept multiple inputs so we can filter on several criteria at once.
    #parser.add_argument("-c", "--column", help="The index of the column to query (starting from 0).")
    #parser.add_argument("-q", "--query", help="The text to match in the specified column.")
    args = parser.parse_args()
    # Validate our input
    if not os.path.isfile(args.input):
        sys.exit("Provided file does not exist: {}".format(args.input))
    main(args)

