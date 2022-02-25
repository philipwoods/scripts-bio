#!/usr/bin/env python
# -*- coding: utf-8
def main(args):
    import pandas as pd

    input_file = args.input_file
    output_file = args.output_file
    data = pd.read_csv(input_file, sep='\t', index_col=False)
    name_column = list(data.columns)[0]
    new_data = data.groupby(name_column).max()
    new_data.to_csv(output_file, sep='\t')

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Fix duplicate functions in ad hoc functional pangenome')
    parser.add_argument('input_file', metavar='FILE', type=str, help='functional occurrence table created with anvi-get-enriched-functions-per-pan-group and --functional-occurrence-table-output option')
    parser.add_argument('output_file', metavar='FILE', type=str, help='functional occurrence table created with anvi-get-enriched-functions-per-pan-group and --functional-occurrence-table-output option')

    args = parser.parse_args()
    main(args)
