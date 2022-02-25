#!/usr/bin/env python3
## Import required modules
import sys
import os
import string

def main():
    ## Global variables
    help_string = """
    Synopsis:
        python evaluate_alignment.py [-h] [--version] [-t <threshold>] [-T] [-nj] input-file

    Description:
        This utility evaluates the proportion of columns in a concatenated alignment
        that are filled for each genome in the alignment. It also recommends genomes
        for exclusion from further analyses using the alignment based on a threshold
        provided by the user. The default threshold is 50%, based on the criteria in
        Parks et al. 2018.

    Options:
        -h | --help
            Display this help text and exit.
        --version
            Display version information and exit.
        -t | --threshold <percentage>
            Specify a minumum percent of columns filled. Default is 50.
        -T | --tab_table
            Output as a tab-separated values table. Reduces memory footprint.
        -nj | --no_judge
            Output won't include the Exclude recommendation column.
    """
    version_string = """
    Last updated July 1 2019
    """

    threshold = 50
    tab_format = False
    judge = True

    ## Handle arguments
    args = sys.argv[1:]
    while len(args) > 0 and args[0][0] == "-":
        arg = args.pop(0)
        if arg in ["-", "--"]:
            break
        elif arg == "-h" or arg == "--help":
            print(help_string)
            sys.exit(0)
        elif arg == "--version":
            print(version_string)
            sys.exit(0)
        elif arg == "-t" or arg == "--threshold":
            threshold = args.pop(0)
            if not is_number(threshold): sys.exit("Must specify a number after " + arg)
            threshold = float(threshold)
        elif arg == "-T" or arg == "--tab_table":
            tab_format = True
        elif arg == "-nj" or arg == "--no_judge":
            judge = False

    if len(args) > 1:
        sys.exit("Extra argument(s) detected.\nUse the --help option to learn more.")

    target = args.pop(0)
    if not os.path.isfile(target):
        sys.exit("Specified file does not exist: " + target)

    parse_alignment(target, threshold, tab_format, judge)

def parse_alignment(target, threshold, tab_format, judge):
    # Initialize variables
    target_file = open(target, 'r')
    columns = ["Genome ID", "Alignment length", "Filled columns", "% filled"]
    if judge: columns.append("Exclude")
    max_id_len = len(columns[0])
    output = []
    output_line = [None, 0, 0, 0]
    if judge: output_line.append(None)
    exclusion_counter = 0

    # Iterate through lines in contcatenated alignment file
    for line in target_file:
        # If the line is a defline...
        if line[0] == ">":
            # Do final calculations, add line to output buffer, and then reset the variable
            if output_line[0] is not None:
                output_line[3] = 100*float(output_line[2])/output_line[1]
                if judge:
                    exclude = meets_exclusion_criteria(output_line[3], threshold)
                    if exclude: exclusion_counter += 1
                    output_line[4] = str(exclude)
                output.append(output_line)
                output_line = [None, 0, 0, 0]
                if judge: output_line.append(None)
            # Split the line at the first space, then remove the first character (">") from the front.
            genome_name = line.split(" ", 1)[0][1:]
            output_line[0] = genome_name
            max_id_len = max(len(genome_name), max_id_len)
        # Otherwise (if it is an alignment line)...
        else:
            # Remove leading and trailing whitespace, including newline characters
            sequence = line.strip()
            # Get the total number of columns
            output_line[1] += len(sequence)
            # Get the number of filled columns
            filled_cols = sequence.replace("-", "")
            output_line[2] += len(filled_cols)
    # Flush out the last output line
    output_line[3] = 100*float(output_line[2])/output_line[1]
    if judge:
        exclude = meets_exclusion_criteria(output_line[3], threshold)
        if exclude: exclusion_counter += 1
        output_line[4] = str(exclude)
    output.append(output_line)
    # Close the concatenated alignment file
    target_file.close()

    # Print output in specified format
    if tab_format:
        print("\t".join(columns))
        for line in output:
            print("\t".join([str(x) for x in line]))
    else:
        # Put our maximum columns widths into a list for convenience
        max_lengths = [len(x) for x in columns]
        max_lengths[0] = max_id_len
        # We want to have three characters of buffer space between each column
        column_spacing = 3
        column_pad = " " * column_spacing
        linewidth = sum(max_lengths) + column_spacing * (len(columns) - 1)
        # We want to left-align the first column, then center the others
        centering = [True]*len(columns)
        centering[0] = False

        # Create pretty table header
        print("-"*linewidth)
        print("Target file: " + os.path.abspath(target))
        if judge:
            print("Min % filled columns: " + str(threshold))
        print("-"*linewidth)
        padded_header = []
        for i in range(len(columns)):
            padded_header.append(pad_text(columns[i], max_lengths[i], centering[i]))
        print(column_pad.join(padded_header))
        print("-"*linewidth)

        # Format output stats
        for line in output:
            padded_line = []
            for i in range(len(columns)):
                text = str(line[i])
                if is_number(text) and not text.isdigit():
                    text = pad_float(text)
                padded_line.append(pad_text(text, max_lengths[i], centering[i]))
            print(column_pad.join(padded_line))

        # Generate summary stats
        print("-"*linewidth)
        print("Total assemblies: " + str(len(output)))
        if judge:
            print("Assemblies to exclude: " + str(exclusion_counter))
            print("Assemblies to keep: " + str(len(output) - exclusion_counter))
        print("-"*linewidth)

def meets_exclusion_criteria(metric, threshold):
    return (metric < threshold)

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
    # Assumes we are dealing with numbers between -99.99 and 999.99
    return "{:6.2f}".format(number)

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
