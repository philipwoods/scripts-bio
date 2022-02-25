#!/usr/bin/env python3
## Import required modules
import os
import sys
import pandas as pd

def main():
    help_string = """
    Synopsis:
        python filter_gc_table.py [-h] [-v] [-l] [-T] [--pretty] [-c COLUMNS] [-r ROWS] GC_TABLE

    Description:
        Takes a positional argument which is a gene cluster presence-absence matrix like
        the output of `anvi-export-table --table gene_cluster_presence_absence PAN_DB`.
        The rows search argument should be a file containing gene cluster names to
        consider when searching. If a rows search file is provided, any entry not listed
        in that file will be excluded from the output. The columns search file should
        have one condition per line specifying allowed column values in the output. Each
        line should follow pandas query() format, e.g. COL1 >= 0. Not all column names
        need to be represented in the file. Even if the input does not have a 'sum'
        column, this can be used as a condition in the columns search file, e.g. sum > 3.

    Options:
        -h | --help
            Display this help text and exit.
        -v | --version
            Display version information and exit.
        -l | --list
            List the column headers and exit.
        -T | --transpose
            Transpose the output.
        -c | --columns
            Provide file containing rules for filtering by column entries.
        -r | --rows
            Provide file containing rows to include in the search.
        --pretty
            Print in a more human-readable format.
    """
    version_string = """
    Last modified 5 September 2019
    """

    list_flag = False
    rows_file = None
    cols_file = None
    pretty_flag = False
    transpose = False

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
        elif arg == "-l" or arg == "--list":
            list_flag = True
        elif arg == "-T" or arg == "--transpose":
            transpose = True
        elif arg == "-c" or arg == "--columns":
            cols_file = args.pop(0)
        elif arg == "-r" or arg == "--rows":
            rows_file = args.pop(0)
        elif arg == "--pretty":
            pretty_flag = True
    if len(args) < 1:
        sys.exit("Missing table input.\nUse the --help option to learn more.")

    gc_table_path = args.pop(0)
    for f in [gc_table_path, rows_file, cols_file]:
        if f != None and not os.path.isfile(f):
            sys.exit("Specified file does not exist: " + f)
    gc_df = pd.read_csv(gc_table_path, sep='\t')

    if list_flag:
        for i in gc_df.columns[1:]:
            print(i)
        sys.exit(1)

    if len(args) > 0:
        sys.exit("Extra argument(s) detected.\nUse the --help option to learn more.")

    out_df = gc_df
    if "sum" not in out_df.columns:
        out_df['sum'] = out_df[out_df.columns[1:]].sum(axis=1)

    # Exclude disallowed rows
    if rows_file != None:
        with open(rows_file, 'r') as f:
            row_list = [x.rstrip() for x in f.readlines()]
            col0 = out_df.columns[0]
            out_df = out_df[out_df[col0].isin(row_list)]

    # Build pandas query string
    if cols_file != None:
        query = []
        with open(cols_file, 'r') as f:
            for line in f:
                query.append("(" + line.rstrip() + ")")
        query_str = " & ".join(query)
        out_df = out_df.query(query_str)

    if transpose:
        out_df = out_df.transpose()

    if pretty_flag:
        print(out_df.to_string(index=transpose))
    else:
        print(out_df.to_csv(sep='\t', index=transpose))

if __name__ == "__main__":
    main()
