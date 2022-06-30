#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt

def main(args):
    # Sequence files are not big enough to be concerned about reading into memory
    with open(args.input, 'r') as f:
        lines = f.readlines()
    # All of the useful information is on the second to last line
    # We split this line by tab, then select the last line to get the actual data
    # We split the data by | to get separate strings for each site in the sequence
    lines = lines[-2].split("\t")[-1].split("|")[:-1]
    # Format data nicely for pandas
    data = defaultdict(list)
    for line in lines:
        fields = line.split(",")
        data['conservation'].append(float(fields[0]))
        data['label'].append(fields[1])
        data['color'].append(fields[-1].strip("[]"))
        if len(fields) == 3:
            data['properties'].append("")
        else:
            data['properties'].append(fields[2].strip())
    # Create dataframe and describe the data
    # If the label is '-' that means the site had more than 25% gaps, so we will exclude those sites
    df = pd.DataFrame(data)
    df_adj = df[df['label'] != "-"]
    print(df_adj.describe().to_string())
    # Plot distributions
    if args.hist is not None:
        fig, ax  = plt.subplots(1)
        bins = [-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5]
        df_adj.hist(column='conservation', bins=bins, ax=ax)
        fig.savefig(args.hist, transparent=False)
    # Write output
    if args.output is not None:
        df_adj.to_csv(args.output, index_label="site")

if __name__ == "__main__":
    desc = ("Reads in a conservation series from Jalview and creates a CSV containing the conservation "
            "scores and other information for each site in the parent alignment.")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('input', help="File containing AA conservation exported in Jalview format.")
    parser.add_argument('-o', '--output', nargs='?', default='output.txt', help="Specify location to write a table of site conservation.")
    parser.add_argument('--hist', default='histogram.png', nargs='?', help="Specify location to write a histogram of conservation scores.")
    args = parser.parse_args()
    if not os.path.isfile(args.input):
        sys.exit("Specified file does not exist: {}".format(args.input))
    main(args)

