#!/usr/bin/env python3
import os
import sys
import argparse
import json
import pandas as pd
from collections import defaultdict
from collections import namedtuple

def main(args)
    return


if __name__ == "__main__":
    desc = ("A command-line utility to view and navigate downloaded BLAST results files. "
            "Your BLAST results must be downloaded in Single-File JSON format from the NCBI results page."
            "Requires pandas to run.")
    epil = ("Note that if a report title contains spaces, you will need to put it in quotes 'like this' when "
            "passing arguments to the --report option.")
    parser = argparse.ArgumentParser(description=desc, epilog=epil)
    parser.add_argument("-i", "--input", required=True, help="The Single-File JSON output from NCBI BLAST.")
    parser.add_argument("-r", "--report", nargs='*', help="Select one or more reports to view in detail. To view all reports, use this option with no arguments.")
    parser.add_argument("--hit", nargs='+', help="Specify the index of one or more hits to display per report. By default, all hits are displayed.")
    parser.add_argument("-a", "--alignments", action='store_true', help="Display alignment reports instead of a summary table.")
    parser.add_argument("-p", "--parameters", action='store_true', help="Display BLAST input parameters in the header of the results.")
    parser.add_argument("-t", "--taxid", action='store_true', help="Display NCBI taxid in the results summary table.")
    parser.add_argument("-l", "--align-length", action='store_true', help="Display alignment length in the results summary table.")
    args = parser.parse_args()
    if not os.path.isfile(args.input):
        sys.exit("Specified file does not exist: {}".format(args.input))
    main(args)

