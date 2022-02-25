#!/usr/bin/env python3
## Import required modules
import os
import sys
import numpy as np
import pandas as pd

def main():
    help_string = """
    This is an accessory script to anvi-script-frequency-enrichment.sh
    The frequency file must be produced from parse-function-frequency.py
    The group file should be formatted like an Anvio layers additional data file with a column labeled 'group' specifying group membership."
    The count file should be formatted with two columns: 'genome' and 'N' where 'N' is the total genes or annotations per genome.
    """
    version_string = """
    Last updated 4 November 2021
    """
    
    frequency_file_path = None
    group_file_path = None
    counts_file_path = None
    
    ## Handle arguments
    args = sys.argv[1:]
    while len(args) > 0 and args[0][0] == "-":
        arg = args.pop(0)
        if arg in ["-", "--"]:
            break
        elif arg == "-h" or arg == "--help":
            print(help_string)
            sys.exit(1)
        elif arg == "-v" or arg == "--version":
            print(version_string)
            sys.exit(1)
        elif arg == "-f" or arg == "--frequency":
            frequency_file_path = args.pop(0)
        elif arg == "-g" or arg == "--groups":
            group_file_path = args.pop(0)
        elif arg == "-c" or arg == "--counts":
            counts_file_path = args.pop(0)

    for f in [frequency_file_path, group_file_path, counts_file_path]:
        if f is None:
            sys.exit("You didn't specify a required file.\nUse the --help option to learn more.")
        elif not os.path.isfile(f):
            sys.exit("Specified file does not exist: " + f)
    if len(args) > 0:
        sys.exit("Extra argument(s) detected.\nUse the --help option to learn more.")

    # Use the group file to create a dictionary:
    # { group1: [sample1, sample2, ...], group2: [sample3, sample4, ...], ...}
    # and column names for the output file.
    groups_df = pd.read_csv(group_file_path, sep='\t')
    group_dict = {}
    for group in groups_df['group'].unique():
        group_dict[group] = groups_df[groups_df['group'] == group]['samples'].to_list()
    
    # Load the frequency file and the counts file, and transpose
    # the counts file so the genome names are column headers.
    frequency_df = pd.read_csv(frequency_file_path, sep='\t')
    counts_df = pd.read_csv(counts_file_path, sep='\t', index_col='genome').T
    # Create an output DataFrame and fill in the columns.
    out_df = pd.DataFrame(columns=['samples'])
    out_df['accession'] = frequency_df['accession']
    out_df['function'] = frequency_df['function']
    # Reverse the column order
    out_df = out_df[out_df.columns[::-1]]
    # Add a column listing the groups associated with each function.
    out_df["associated_groups"] = "placeholder"
    # Add columns for the sample size and proportion of each group with the function.
    for group,members in group_dict.items():
        # The .sum() call returns a pandas Series, but in the case of the counts df,
        # we know this will always have only a single value because it only has one row.
        # We need to extract the value before dividing because Series / Series doesn't
        # work when they have different lengths.
        grp_count = counts_df[members].sum(axis=1)[0]
        out_df["p_{}".format(group)] = frequency_df[members].sum(axis=1) / grp_count
        out_df["N_{}".format(group)] = grp_count
    # Use the proportions and counts just added to fill the associated_groups column.
    for accession in out_df['accession']:
        # Pandas treats entries in a list as the values for each row in a column. This
        # means that casting this df to an array gives [[a, b, c, d, ...]] since we are
        # only pulling one row (one accession).
        p_vector = np.array(out_df[out_df['accession'] == accession].filter(regex='p_'))[0]
        N_vector = np.array(out_df[out_df['accession'] == accession].filter(regex='N_'))[0]
        # For each accession, this takes in two vectors
        #           grp1,   grp2,   grp3, ...
        # p_vect = [ p1,     p2,     p3,  ... ]
        # N_vect = [ N1,     N2,     N3,  ... ]
        # And outputs one vector
        # output = [ b1,     b2,     b3,  ... ]
        # where the booleans are True or False depending on whether the corresponding group is enriched
        enriched_groups_vector = get_enriched_groups(p_vector, N_vector)
        # We get the group labels of the enriched groups
        associated_groups = [c for i,c in enumerate(group_dict.keys()) if enriched_groups_vector[i]]
        if not associated_groups: # If there aren't any associated groups...
            associated_groups.append("NA")
        # Then we insert the proper value in the proper position.
        # See comment on line 77 for why the pandas query looks like this.
        function_index = out_df.index[out_df['accession'] == accession].to_list()[0]
        out_df.loc[function_index,'associated_groups'] = ",".join(associated_groups)
    # Write to output file.
    print(out_df.to_csv(sep='\t', index=False))

# Taken from anvio utils.py
def get_enriched_groups(props, reps):
    '''
        Accepts a vector of proportions and number of replicates per group and
        returns a boolean vector where each group that has proportion above
        the "expected" (i.e. the overall proportion) is True and the rest are False.
    '''
    # if the function doesn't occur at all then test_statistic is zero and p-value is 1
    if not np.count_nonzero(props):
        return np.zeros(len(props))
    overall_portion = np.sum(np.multiply(props, reps)) / np.sum(reps)

    return props > overall_portion

if __name__ == "__main__":
    main()

