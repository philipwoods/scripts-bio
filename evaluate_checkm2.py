#!/usr/bin/env python3
## Import required modules
import sys
import os
import pandas as pd
import argparse

def main(args):
    df = pd.read_table(args.report, header=0, index_col=False)
    df["Quality"] = df["Completeness"] - 5 * df["Contamination"]
    df["Passed"] = df["Quality"] > args.min_quality
    df.to_csv(args.output, sep='\t', columns=["Name", "Completeness", "Contamination", "Quality", "Passed"], index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', '-o', default=sys.stdout, help="Location to write output. Default: stdout")
    parser.add_argument('--min-quality', default=50, type=float, help="Threshold for minimum acceptable quality score. Default: 50")
    parser.add_argument('report', help="Quality report output from CheckM2.")
    args = parser.parse_args()
    if not os.path.isfile(args.report):
        sys.exit("The provided file does not exist: {}".format(args.quality_report))
    main(args)
