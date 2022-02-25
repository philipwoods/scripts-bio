#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd

def main():
    # Handle arguments
    function_desc = ("Condense multiple one-dimensional data files with a common index into a two-dimensional file. "
                     "Each index value does not need to be present in all input files. "
                     "Uses input file names to label entries in the output.")
    parser = argparse.ArgumentParser(allow_abbrev=False, description=function_desc)
    parser.add_argument("-c", "--cumulative", action='store_true',
                        help="Report cumulative data, i.e. reverse cumulative sum along the index for each file. Default is off.")
    parser.add_argument("--no-header", dest='header', default=0, action='store_const', const=None,
                        help="Indicate that the input files do not have a header row. By default, assumes the first row is a header.")
    parser.add_argument("--fill-index", action='store_true',
                        help="Indicate that the index is of integer type, and that any missing intermediate indices should be filled in in the output.")
    parser.add_argument("--fill-data", metavar="FILL",
                        help="The value to be used to fill missing data values. Left empty by default.")
    parser.add_argument("--head", help="Substring to remove from the start of each filename when creating labels for the final data.")
    parser.add_argument("--tail", help="Substring to remove from the end of each filename when creating labels for the final data.")
    parser.add_argument("-i", "--index", default="index",
                        help="A name to be used in the output describing the index column in the input files. Default is 'index'.")
    parser.add_argument("-s", "--separator", default='\t',
                        help="The character separating the columns in the input files. Default is TAB.")
    parser.add_argument("FILE", nargs="+", help="One-dimensional data file containing an index column followed by a numerical data column.")
    args = parser.parse_args()

    # Check validity of arguments
    files = [f for f in args.FILE if os.path.isfile(f)]
    if len(args.FILE) != len(files):
        notfiles = [os.path.abspath(f) for f in args.FILE if f not in files]
        sys.exit("At least one provided file does not exist:\n{}".format("\n".join(notfiles)))
    if args.cumulative and not args.fill_data:
        print("WARNING: Using the --cumulative flag without the --fill-data option may create missing data in your output.")

    # Begin importing input files
    data = []
    for f in files:
        filename = os.path.basename(f)
        if args.head and filename.startswith(args.head):
            size = len(args.head)
            filename = filename[size:]
        if args.tail and filename.endswith(args.tail):
            size = len(args.tail) * -1
            filename = filename[:size]
        series = pd.read_csv(f, sep=args.separator, header=args.header, names=[args.index, filename], index_col=args.index)
        data.append(series)
    # Condense the input data
    out_df = pd.concat(data, axis=1, sort=True).convert_dtypes()
    if args.fill_index:
        max_index = max(out_df.index)
        min_index = min(out_df.index)
        out_df = out_df.reindex(range(min_index, max_index + 1)).convert_dtypes()
    if args.fill_data:
        out_df = out_df.fillna(pd.to_numeric(args.fill_data)).convert_dtypes()
    if args.cumulative:
        # Pandas cumulative sum goes in the reverse direction from what we want.
        # We solve this by reversing the index order before the sum, then reversing again afterwards.
        acc_df = out_df.loc[::-1].cumsum()[::-1]
        out_df = acc_df
    # Transpose the DataFrame before output
    print(out_df.T.to_csv(index_label=args.index, sep='\t'))

if __name__ == "__main__":
    main()

