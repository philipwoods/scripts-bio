#!/usr/bin/env python3
### Import required modules
import os
import sys
import pandas as pd

# TO DO: the ngram step takes really long, is there a more efficient way to do that?

def main():
    help_string = """
    Synopsis:
        python deduplicate-synteny.py SYNTENY_FILE OUT_FILE
    
    Description:
        This is an utility for the script anvi-script-extract-synteny-groups.sh
        and is not generally intended for independent use. For that script to work,
        both scripts must be in the same directory. Assumes the input file has the
        following format:

        ngram    N   genome   count

        Each line holds a record of a set of ordered genes in a single genome.
        The ngram is the description of the genes in order.
        N is the length of the ngram (how many genes in order).
        The genome is the genome this ngram occurs in.
        The count is the number of times it occurs in the specified genome.
    """
    version_string = """
    Last updated 12 November 2021
    """
    ## Handle arguments
    if '-h' in sys.argv or '--help' in sys.argv:
        print(help_string)
        sys.exit(1)
    if '-v' in sys.argv or '--version' in sys.argv:
        print(version_string)
        sys.exit(1)
    synteny_path = sys.argv[1]
    out_path = sys.argv[2]

    # Read in the data.
    synteny_df = pd.read_csv(synteny_path, names=['ngram','N','genome','count'], sep='\t', dtype={'N':"Int8", 'count':"Int8"}, skip_blank_lines=False)
    genomes = synteny_df['genome'].dropna().unique()
    for genome in genomes: # For each genome...
        print("Analyzing {}...".format(genome))
        ngrams = synteny_df[synteny_df['genome'] == genome]['ngram'].dropna().unique()
        for ngram in ngrams: # and each ngram...
            # Get all rows which:
            # 1) are from the right genome,
            # 2) are not the ngram of interest,
            # 3) have an ngram that contains the ngram of interest.
            # Then use the 'N' column to see if there are any such rows.
            sub = synteny_df[(synteny_df['genome'] == genome) & ~(synteny_df['ngram'] == ngram) & (synteny_df['ngram'].str.contains(ngram, regex=False))]['N'].sum()
            if sub > 0: # If there are any such rows, decrement the count for that ngram in that genome.
                synteny_df.loc[((synteny_df['ngram'] == ngram) & (synteny_df['genome'] == genome)), 'count'] -= 1
    # Once this process is finished, remove any ngrams without a positive count.
    duplicates = synteny_df[synteny_df['count'] < 1].index
    synteny_df.drop(index=duplicates, inplace=True)
    synteny_df.to_csv(out_path, sep='\t', index=False, header=False)
    
if __name__ == "__main__":
    main()

