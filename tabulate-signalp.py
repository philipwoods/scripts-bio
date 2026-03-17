#!/usr/bin/env python3
import os
import sys
import argparse
import json
import pandas as pd
from collections import defaultdict
from collections import namedtuple

def main(args):
    # Load JSON input
    with open(args.input, 'r') as f:
        results = json.load(f)
    data_dict = results['SEQUENCES']
    out_dict = []
    for entry in data:
        seqname = entry['Name']
        prediction = entry['Prediction']
        prot_types = entry['Protein_types']
        likelihoods = entry['Likelihood']
        row = {
            'seq_id': seqname,
            'prediction': prediction
        }
        for name, p in zip(prot_types, likelihoods):
            row[name] = p
        output.append(row)
    out_df = pd.DataFrame.from_records(out_dict)
    out_df.to_csv(args.output, index=False, sep=args.sep)

if __name__ == "__main__":
    desc = ("A command-line utility to convert the standard output from SignalP from JSON to "
            "a tabular format.")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("input", help="The JSON output from SignalP 6.0")
    parser.add_argument("--sep", default='\t', help="The separator to use in the output table. Default: TAB")
    parser.add_argument("--output", default=sys.stdout, help="File name to write output. Default: STDOUT")
    args = parser.parse_args()
    if not os.path.isfile(args.input):
        sys.exit("Specified file does not exist: {}".format(args.input))
    if len(args.sep) > 1:
        sys.exit("Separator must be exactly one character long.")
    main(args)

