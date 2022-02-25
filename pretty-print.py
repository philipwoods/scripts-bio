#!/usr/bin/env python3
## Import required modules
import sys
import os
import csv

def main():
    ## Global variables
    help_string = """
    Usage:
        python pretty-print.py [-h] [--version] [--no-header] <input-file>
    
    Description:
        Takes in a TAB-delimited file and prints a pretty version formatted
        for easy reading. Assumes the first line is column headers.

    Options:
        -h | --help
            Display this help text and exit.
        --version
            Display version information and exit.
        --no-header
            This flag indicates there are no column headers.
    """
    version_string = """
    Last updated 25 November 2021
    """

    header = True

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
        elif arg == "--no-header":
            header = False
        elif arg not in ["-", "--"]:
            sys.exit("Invalid option: " + arg)

    if not args:
        sys.exit("You must provide an input file.\nUse the --help option to learn more.")
    if len(args) > 1:
        sys.exit("Extra argument(s) detected.\nUse the --help option to learn more.")

    target = args.pop(0)
    if not os.path.isfile(target):
        sys.exit("Specified file does not exist: " + target)

    ## Parse and print the file
    with open(target,'r') as f:
        filereader = csv.reader(f, delimiter='\t', skipinitialspace=True)
        filelist = list(filereader)

        # We need to know the longest entry in each column
        ncols = len(filelist[0])
        maxlengths = []
        for i in range(ncols):
            # Need to discard trailing empty lines
            maxlengths.append(max([len(line[i]) for line in filelist if line]))
        # We want to have three characters of buffer space between each column
        column_spacing = 3
        column_pad = " " * column_spacing
        linewidth = sum(maxlengths) + column_spacing * (ncols - 1)
        # We want to center all of the columns except the first
        centering = [True] * ncols
        centering[0] = False

        # Generate lines into an output buffer
        output = []
        for line in filelist:
            paddedline = []
            for i in range(ncols):
                text = ""
                if line: text = line[i]
                if not text.isdigit() and is_number(text):
                    text = pad_float(text)
                paddedline.append(pad_text(text, maxlengths[i], centering[i]))
            output.append(column_pad.join(paddedline))
        if header:
            output.insert(0, "-"*linewidth)
            output.insert(2, "-"*linewidth)
        # Print the buffer
        print("\n".join(output))

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

if __name__ == "__main__":
    main()

