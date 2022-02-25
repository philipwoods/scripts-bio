#!/usr/bin/env python3
import os
import sys
import argparse
import tempfile
import pandas as pd

def main(args):
    # Declare a header and initialize the output buffer
    header = ['gene_cluster_id']
    header.extend(args.annotation_sources)
    out_buffer = ['\t'.join(header)]
    # Add information for each requested cluster
    for i, cluster in enumerate(args.clusters):
        print("Working on cluster {c} ({n} of {t})".format(c=cluster, n=i, t=len(args.clusters)))
        out_line = [cluster]
        temp = os.path.join(tempfile.gettempdir(), "temp_{}".format(cluster))
        # Grep for the cluster to reduce memory cost
        os.system("head -1 {data} > {tmpfile}".format(data=args.summary, tmpfile=temp))
        os.system("grep -w {pattern} {data} >> {tmpfile}".format(pattern=cluster, data=args.summary, tmpfile=temp))
        filtered_summary = pd.read_csv(temp, sep='\t', index_col='unique_id')
        total = filtered_summary['gene_cluster_id'].count() # Total number of times the cluster appears
        filtered_summary = filtered_summary.fillna("NO ANNOTATION")
        # Add a column for each annotation source
        for source in args.annotation_sources:
            sizes = filtered_summary[source].groupby(filtered_summary[source]).size().sort_values(ascending=False)
            functions = ["{a} (N={n}/{t})".format(a=annotation, n=num, t=total) for annotation, num in zip(sizes.index, sizes.values)]
            out_line.append(", ".join(functions))
        # Append this cluster's info to the output buffer
        out_buffer.append('\t'.join(out_line))
        os.remove(temp)
    # Write output
    with open(args.output, 'w') as f:
        f.write('\n'.join(out_buffer))

if __name__ == "__main__":
    # Handle inputs
    desc = ("Provides information about the functional annotations of one or more gene clusters in a pangenome. "
            "Be careful: Depending on the settings used to create the pangenome, it may or may not be correct to "
            "assume that all genes within a gene cluster actually have the same function.")
    epil = ("The output format will have one line per gene cluster and one column per annotation source. "
            "The entries for each annotation source will be formatted as 'ANNOTATION (N=ANNOTATION_FREQ/CLUSTER_FREQ),...' "
            "where ANNOTATION_FREQ is the number of times that gene cluster has the given annotation and CLUSTER_FREQ is "
            "the number of genes in the pangenome which are part of that cluster.")
    parser = argparse.ArgumentParser(description=desc, epilog=epil)
    parser.add_argument("-s", "--summary", required=True,
        help="The output of anvi-summarize, describing gene clusters in a pangenome.")
    parser.add_argument("-c", "--clusters", required=True, nargs='+', metavar="CLUSTER",
        help="The ID of the gene cluster(s) to be summarized.")
    parser.add_argument("-a", "--annotation-sources", required=True, nargs='+', metavar="SOURCE",
        help="The annotation source(s) to report summary information from.")
    parser.add_argument("-o", "--output", required=True, help="Path for output to be written.")
    args = parser.parse_args()
    # Check for errors
    if not os.path.isfile(args.summary):
        sys.exit("File does not exist: {}".format(args.summary))
    # Run main function
    main(args)

