#!/usr/bin/env python3
## Import required modules
import sys
import os
import csv

def main():
    ## Global variables
    help_string = """
    Synopsis:
        python evaluate_qa.py [-h] [--version] [-T] [<threshold options>] input-file

    Description:
        This utility helps extract quality criteria from an output file in tab table
        format produced by the checkm qa command using the --tab_table option. It
        considers a quality metric based on completeness and contamination, the
        number of scaffolds in an assembly, and the N50 of those scaffolds. If any
        of these criteria do not meet the required threshold, it is recommended that
        the assembly not be used in further analysis.

        The quality metric is (% completeness - 5 * % contamination). This metric and
        the defaults for each criterion are based on Parks et al. (2017) Nature
        Microbiology 2, 1533--1542.
    
    Options:
        -h | --help
            Display this help text and exit.
        --version
            Display version information and exit.
        -T | --tab_table
            Output as a tab-separated values table. Reduces memory footprint.
        -qt | --quality_threshold <min-quality>
            Specify a minumum quality threshold. Default is 50.
        -ms | --max_scaffolds <max-scaffolds>
            Specify a maximum number of scaffolds. Default is 500.
        -n50 | --min_n50 <min-N50>
            Specify a minumum N50 for the scaffolds. Default is 10kb.
    """
    version_string = """
    Last updated June 20 2019
    """

    tab_format = False
    thresholds = {
        'min_quality':      50,
        'max_scaffolds':    500,
        'min_n50':          10000
    }

    ## Handle arguments
    args = sys.argv[1:]
    while len(args) > 0 and args[0][0] == "-":
        arg = args.pop(0)
        if arg == "-h" or arg == "--help":
            print(help_string)
            sys.exit(0)
        elif arg == "--version":
            print(version_string)
            sys.exit(0)
        elif arg == "-T" or arg == "--tab_table":
            tab_format = True
        elif arg == "-qt" or arg == "--quality_threshold":
            thresholds['min_quality'] = args.pop(0)
            if not is_number(thresholds['min_quality']): sys.exit("Must specify a number after " + arg)
            thresholds['min_quality'] = float(thresholds['min_quality'])
        elif arg == "-ms" or arg == "--max_scaffolds":
            thresholds['max_scaffolds'] = args.pop(0)
            if not is_number(thresholds['max_scaffolds']): sys.exit("Must specify a number after " + arg)
            thresholds['max_scaffolds'] = float(thresholds['max_scaffolds'])
        elif arg == "-n50" or arg == "--min_n50":
            thresholds['min_n50'] = args.pop(0)
            if not is_number(thresholds['min_n50']): sys.exit("Must specify a number after " + arg)
            thresholds['min_n50'] = float(thresholds['min_n50'])
        elif arg not in ["-", "--"]:
            sys.exit("Invalid option: " + arg)

    if len(args) > 1:
        sys.exit("Extra argument(s) detected.\nUse the --help option to learn more.")

    target = args.pop(0)
    if not os.path.isfile(target):
        sys.exit("Specified file does not exist: " + target)

    parse_qa(target, tab_format, thresholds)

def parse_qa(target_file, tab_format, thresholds_dict):
    with open(target_file,'r') as f:
        qa_reader = csv.reader(f, delimiter='\t', skipinitialspace=True)

        # Identify the columns of interest
        header = next(qa_reader)
        # Columns are listed in the order we want to output them.
        columns = ["Bin Id", "Marker lineage", "# scaffolds", "N50 (scaffolds)", "Completeness", "Contamination"]
        col_indices = [header.index(x) for x in columns]
        # Add headers for columns that we'll calculate from the file data
        columns.append("Quality")
        columns.append("Exclude")

        # Construct output according to flagged format
        if tab_format:
            print('\t'.join(columns))
            for line in qa_reader:
                # Calculate metrics and compare to thresholds
                complete_i = col_indices[columns.index("Completeness")]
                contam_i = col_indices[columns.index("Contamination")]
                scaffolds_i = col_indices[columns.index("# scaffolds")]
                n50_i = col_indices[columns.index("N50 (scaffolds)")]
                metrics = {
                    'quality': quality_metric(line[complete_i], line[contam_i]),
                    'scaffolds': float(line[scaffolds_i]),
                    'n50': float(line[n50_i])
                }
                exclude = meets_exclusion_criteria(metrics, thresholds_dict)
                # Print output
                output = [line[i] for i in col_indices]
                output.append(str(metrics['quality']))
                output.append(str(exclude))
                print('\t'.join(output))
        else:
            # Initialize an output string 'buffer' where each element is one line of output
            output = []

            # Find the max length of strings in different columns for pretty formatting.
            qa_list = list(qa_reader)  # Convert reader to a list of file lines.
            qa_list.insert(0, header)  # Add the header back to get the complete file.
            max_lengths = []
            # We look through the file to find the longest string length in each
            # of our columns of interest.
            for i in col_indices:
                max_lengths.append(max([len(line[i]) for line in qa_list]))
            # For the added columns, the entries are at most 6-character numbers
            # so the column headers will be the longest strings.
            for col_name in columns[len(col_indices):]:
                max_lengths.append(len(col_name))

            # We want to have three characters of buffer space between each column
            column_spacing = 3
            column_pad = " " * column_spacing
            linewidth = sum(max_lengths) + column_spacing * (len(columns) - 1)

            # Create pretty table header
            output.append("-"*linewidth)
            output.append("Target file: " + os.path.abspath(target_file))
            output.append("Min quality: " + str(thresholds_dict['min_quality']))
            output.append("Max # scaffolds: " + str(thresholds_dict['max_scaffolds']))
            output.append("Min scaffold N50: " + str(thresholds_dict['min_n50']) + "bp")
            centering = [True]*len(columns) # Set a centering mask for the columns
            centering[0] = False    # First column (Bin Id) should be left justified
            output.append("-"*linewidth)
            padded_header = []
            for i in range(len(columns)):
                padded_header.append(pad_text(columns[i], max_lengths[i], centering[i]))
            output.append(column_pad.join(padded_header))
            output.append("-"*linewidth)

            # Add the data
            exclusion_counter = 0
            for line in qa_list[1:]:
                # Calculate metrics and compare to thresholds
                complete_i = col_indices[columns.index("Completeness")]
                contam_i = col_indices[columns.index("Contamination")]
                scaffolds_i = col_indices[columns.index("# scaffolds")]
                n50_i = col_indices[columns.index("N50 (scaffolds)")]
                metrics = {
                    'quality': quality_metric(line[complete_i], line[contam_i]),
                    'scaffolds': float(line[scaffolds_i]),
                    'n50': float(line[n50_i])
                }
                exclude = meets_exclusion_criteria(metrics, thresholds_dict)
                if exclude: exclusion_counter += 1

                # Generate line in buffer
                padded_line = []
                for i in range(len(columns)):
                    if i < len(col_indices):
                        text = line[col_indices[i]]
                    elif i == columns.index("Quality"):
                        text = pad_float(metrics['quality'])
                    elif i == columns.index("Exclude"):
                        text = str(exclude)
                    if not text.isdigit() and is_number(text):
                        text = pad_float(text)
                    padded_line.append(pad_text(text, max_lengths[i], centering[i]))
                # Add the finished line to the buffer
                output.append(column_pad.join(padded_line))

            # Generate summary stats
            output.append("-"*linewidth)
            output.append("Total assemblies: " + str(len(qa_list[1:])))
            output.append("Assemblies to exclude: " + str(exclusion_counter))
            output.append("Assemblies to keep: " + str(len(qa_list[1:]) - exclusion_counter))
            output.append("-"*linewidth)
            # Print the output buffer
            print("\n".join(output))

def quality_metric(completeness, contamination):
    # The items provided by the csv reader are strings
    return float(completeness) - 5 * float(contamination)

def meets_exclusion_criteria(metrics, thresholds):
    # Compare each metric to its threshold to see if it fails the criterion
    quality = (metrics['quality'] < thresholds['min_quality'])
    scaffolds = (metrics['scaffolds'] > thresholds['max_scaffolds'])
    n50 = (metrics['n50'] < thresholds['min_n50'])
    # If any of these criteria are failed (True) then it should be excluded.
    return (quality or scaffolds or n50)

def pad_text(text,length,centered):
    pad_length = length - len(text)
    if centered:
        start_pad = " " * (pad_length//2)
        end_pad = " " * (pad_length - pad_length//2)
    else:
        start_pad = ""
        end_pad = " " * pad_length
    return start_pad + text + end_pad

def pad_float(number):
    number = float(number)
    # Format the number as a 7-character-wide string with two decimal places
    # Assumes we are dealing with numbers between -999.99 and 9999.99
    return "{:7.2f}".format(number)

def is_number(s):
    # We have to do this because in python 2 there is no str.isnumeric()
    # As above, makes assumptions about possible inputs e.g. no scientific notation.

    # Format of /^-?\d*\.?\d*$/ with at least one digit somewhere
    condition1 = s.lstrip("-").replace(".","",1).isdigit()
    # No more than one leading '-' character
    condition2 = ((len(s.lstrip("-")) + 1) >= len(s))
    return (condition1 and condition2)

if __name__ == "__main__":
    main()
