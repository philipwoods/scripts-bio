#!/usr/bin/env python3
# Import required modules
import os
import sys
import pandas as pd

def main():
    help_string = """
    Synopsis:
        python parse-function-frequency.py <input-dir> <mode> <frequency-out> <counts-out>
    
    Description:
        This is an auxiliary utility for the script anvi-script-frequency-enrichment.sh
        and is not generally intended for independent use. For example, it has no input
        validation and relies on its parent script always passing well-formed input. For
        that script to work, both scripts must be in the same directory.
        
        The input directory should contain one file per genome listing annotation information.
        The mode parameter allows you to choose how the annotations should be interpreted.
        The available modes are listed below.

        'annotation'    Provide the number of times an annotation is present in each genome.
        'gene'          Provide the number of genes with a given annotation in each genome.

        The last two arguments are the file paths to use for output..
    """
    version_string = """
    Last updated 4 November 2021
    """

    if "-h" in sys.argv or "--help" in sys.argv:
        print(help_string)
        sys.exit(1)

    #Handle arguments
    input_dir = sys.argv[1]
    mode = sys.argv[2]
    freq_out = sys.argv[3]
    count_out = sys.argv[4]

    # Import the data files we need and extract the relevant genome names
    file_list = os.listdir(input_dir)
    file_list = list(filter(lambda x: "-functions.tmp" in x, file_list))
    genomes = [name.replace("-functions.tmp","") for name in file_list]

    # Each input file will have the following format, where each line represents a single gene:
    # accession     function
    # X             textX
    # Y!!!Z         textY!!!textZ

    # The output file will have the following format:
    # accession     function    genome1 genome2 genome3 ...
    # X             textX       N1_X    N2_X    N3_X    ...
    # Y             textY       N1_Y    N2_Y    N3_Y    ...
    # Z             textZ       N1_Z    N2_Z    N3_Z    ...

    # Create a DataFrame to store the frequency output
    cols = ['accession', 'function']
    cols.extend(genomes)
    freq_df = pd.DataFrame(columns=cols)

    # For each genome's file...
    for func_file in file_list:
        filename = input_dir + '/' + func_file
        genome = func_file.replace("-functions.tmp","")
        func_df = pd.read_csv(filename, sep='\t')
        # For each annotated gene...
        for index, row in func_df.iterrows():
            accessions = row['accession'].split("!!!")
            functions = row['function'].split("!!!")
            # Add a row to the output for any annotations not already present
            for accession, function in zip(accessions, functions):
                if accession not in freq_df['accession'].values:
                    newrow = [accession, function]
                    newrow.extend([0]*len(genomes))
                    freq_df = freq_df.append(pd.DataFrame(columns=cols, data=[newrow]))
            # We have to handle annotations differently depending on the output mode
            if mode == 'annotation':
                # Increment the annotation count(s) for this genome by stepping through
                # each annotation for the gene. If an accession is present multiple times,
                # its counter will increment multiple times.
                for accession in accessions:
                    freq_df.loc[freq_df['accession'] == accession, genome] += 1
            elif mode == 'gene':
                # Increment the annotation count(s) for this genome by testing if an
                # accession is present on a gene. If an accession is present multiple times,
                # its counter will increment one time.
                freq_df.loc[freq_df['accession'].isin(accessions), genome] += 1

    # Output final counts of annotations to be used by later scripts.
    # The values of N are only meaningful in 'annotation' mode.
    count_df = pd.DataFrame(columns=['genome', 'N'])
    for genome in genomes:
        count = freq_df[genome].sum()
        count_df = count_df.append(pd.DataFrame([[genome, count]], columns=['genome', 'N']))

    freq_df.to_csv(freq_out, sep='\t', index=False)
    count_df.to_csv(count_out, sep='\t', index=False)

if __name__ == "__main__":
    main()
