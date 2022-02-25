#!/usr/bin/env python3
# Import required modules
import os
import sys
import pandas as pd

def main():
    help_string = """
    Synopsis:
        python parse-cog-cat-functions.py INPUT_DIRECTORY PRETTY_FLAG
    
    Description:
        This is an auxiliary utility for the script anvi-get-cog-category-frequencies.sh
        and is not generally intended for independent use. For example, it has no input
        validation and relies on its parent script always passing well-formed input. For
        that script to work, both scripts must be in the same directory. The pretty flag
        should be passed as 1 for true or 0 for false.
    """
    version_string = """
    Last updated 5 November 2019
    """

    if "-h" in sys.argv or "--help" in sys.argv:
        print(help_string)
        sys.exit(1)

    #Handle arguments
    input_dir = sys.argv[1]
    pretty_flag = bool(int(sys.argv[2]))

    # Import the data files
    file_list = os.listdir(input_dir)
    file_list = list(filter(lambda x: "_functions_COG_CATEGORY.tab" in x, file_list))

    # Initialize useful variables
    categories = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    frequencies = []

    for cog_file in file_list:
        filename = input_dir + '/' + cog_file
        with open(filename, 'r') as f:
            # First line is a column header, then every line has a single letter category.
            cogs = f.readlines()[1:]
        # Remove trailing newlines from each entry
        cogs = [x.strip() for x in cogs]

        genome_name = cog_file.replace("_functions_COG_CATEGORY.tab", "")
        row = [genome_name]
        for category in categories:
            row.append(cogs.count(category))
        row.append(len(cogs))
        frequencies.append(row)

    # Create a DataFrame
    columns = ['genome_name']
    columns.extend(list(categories))
    columns.append("Total")
    df = pd.DataFrame(frequencies, columns=columns)
    if pretty_flag:
        print(df.to_string(index=False))
    else:
        print(df.to_csv(sep='\t', index=False))

if __name__ == "__main__":
    main()
