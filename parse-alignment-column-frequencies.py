#!/usr/bin/env python3
# Import required modules
import os
import sys
from Bio import AlignIO
import pandas as pd

def main():
    help_string = """
    Synopsis:
        python parse-alignment-column-frequencies.py <output-dir> <alignment> <groups-file>
    
    Description:
        This is an auxiliary utility for the script anvi-script-alignment-enrichment.sh
        and is not generally intended for independent use. For example, it has no input
        validation and relies on its parent script always passing well-formed input. For
        that script to work, both scripts must be in the same directory.
    """
    version_string = """
    Last updated 4 November 2021
    """

    if "-h" in sys.argv or "--help" in sys.argv:
        print(help_string)
        sys.exit(1)

    #Handle arguments
    output_dir = sys.argv[1]
    alignment_file= sys.argv[2]
    groups_file = sys.argv[3]

if __name__ == "__main__":
    main()
