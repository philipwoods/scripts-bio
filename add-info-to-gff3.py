#!/usr/bin/env python3
## Import required modules
import os
import sys
import pandas as pd

def main():
    help_string = """
    Synopsis:
        python add-info-to-gff3 <GFF3-feature-file> <tRNAs-file> <gene-cluster-file>
    
    Description:
        This is an auxiliary utility for the script anvi-script-get-annotated-genome.sh
        and is not generally intended for independent use. For example, it doesn't
        have any input validation to speak of, relying on the fact that the script
        above will only pass well-formed input. For that script to work, both scripts
        must be in the same directory.
    """
    version_string = """
    Last updated 23 October 2021
    """

    ## Handle arguments
    if '-h' in sys.argv or '--help' in sys.argv:
        print(help_string)
        sys.exit(1)
    if '-v' in sys.argv or '--version' in sys.argv:
        print(version_string)
        sys.exit(1)
    gff3_path = sys.argv[1]
    trna_path = sys.argv[2]
    gc_path = None
    if len(sys.argv) > 3:
        gc_path = sys.argv[3]
    temp_path = gff3_path + ".tmp"

    # Import the additional data files into a pandas dataframe
    print("Reading in additional data files...")
    trna_df = pd.read_table(trna_path, dtype=str)
    gc_df = None
    if gc_path != None:
        gc_df = pd.read_csv(gc_path, sep='\t', dtype=str)

    # Open the input file and output temp file, then write the header to the output.
    print("Beginning writing process...")
    gff3_in = open(gff3_path, 'r')
    gff3_out = open(temp_path, 'w')
    gff3_out.write(gff3_in.readline())

    # Find additional data for each gene and write the updated line to the temp file.
    for line in gff3_in:
        # As of anvio v7.1, the GFF3 file uses ID=<project-name>___<gene_caller_id> where
        # <project-name> was specified during the creation of the contigs db we are using.
        # In genes with no function, the line ends after the ID. If there is a function,
        # there is a semicolon after the ID.
        start = line.find("___") + 3 # Get the index of the start of the ID number
        end = line.find(";") # Try to get the index of the end of the ID number
        if end > 0: # If there was a semicolon...
            gene_id = line[start:end]
        else: # If there wasn't...
            gene_id = line[start:].strip() # Strip the \n off the end of the line
        # The Transfer_RNAs annotation source includes semicolons in the function, but GFF3 reserves them as special
        trna = "!!!".join(trna_df[trna_df['gene_callers_id']==gene_id]['function'].to_numpy()).replace(";",",")
        accession = "!!!".join(trna_df[trna_df['gene_callers_id']==gene_id]['accession'].to_numpy())
        new_attributes = ""
        if accession: # If there is a value, add fields for it.
            new_attributes = ";Name=" + accession + ";product=" + trna
        # Start to build the output line
        out_line = line.rstrip()
        out_line = out_line + new_attributes
        if gc_path != None:
            gene_cluster = "!!!".join(gc_df[gc_df['gene_caller_id']==gene_id]['gene_cluster_id'].to_numpy())
            out_line = out_line + ";Alias=" + gene_cluster
        out_line = out_line + "\n"
        gff3_out.write(out_line)

    # Close relevant files and replace the old feature annotation file with the new one.
    print("Cleaning up...")
    gff3_in.close()
    gff3_out.close()
    os.replace(temp_path, gff3_path)

if __name__ == "__main__":
    main()
