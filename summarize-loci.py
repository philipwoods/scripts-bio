#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd

def main(args):
    # Convert inputs to useful forms.
    summary = pd.read_csv(args.summary, sep="\t", index_col='unique_id')
    columns = ['gene_cluster_id', 'bin_name', 'genome_name']
    if args.columns:
        columns.extend(args.columns)
    # Handle the two methods of passing locus definitions
    loci = []
    if args.loci:
        with open(args.loci, 'r') as f:
            loci = [line.strip() for line in f]
    else:
        loci.append(args.locus)
    # Construct locus summaries.
    for locus_id, locus in enumerate(loci):
        clusters = locus.split("::")
        # Filter the input summary file.
        out_df = summary.loc[summary['gene_cluster_id'].isin(clusters), columns]
        # Write the output.
        out_path = os.path.join(args.out_dir, "locus-{}.csv".format(locus_id))
        out_df.to_csv(out_path, sep='\t')

if __name__ == "__main__":
    # Handle inputs
    parser = argparse.ArgumentParser(description="Provides information about one or more genomic loci based on the pangenome gene clusters contained within.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--locus", help="Mutually exclusive with --loci. A description of the locus in the format GC_####::GC_#####::...")
    group.add_argument("--loci", help="Mutually exclusive with --locus. A file with each line containing one locus description in the format above.")
    parser.add_argument("-s", "--summary", required=True, help="The output of anvi-summarize, describing gene clusters in a pangenome.")
    parser.add_argument("-c", "--columns", nargs='+',
                        help="Columns from the summary file to report in the output. The columns gene_cluster_id, bin_name, and genome_name are always reported.")
    parser.add_argument("-o", "--out-dir", required=True, help="Directory path for output to be written. One file will be created for each locus.")
    args = parser.parse_args()
    # Check for errors
    if not os.path.isfile(args.summary):
        sys.exit("File does not exist: {}".format(args.summary))
    if args.loci and not os.path.isfile(args.loci):
        sys.exit("File does not exist: {}".format(args.loci))
    # Run main function
    main(args)

