#!/usr/bin/env python3
import os
import sys
import argparse
import re
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

def make_lalign_formatter(df, cols=None):
    """
    Construct formatter dict to left-align columns.

    Parameters
    ----------
    df : pandas.core.frame.DataFrame
        The DataFrame to format
    cols : None or iterable of strings, optional
        The columns of df to left-align. The default, cols=None, will
        left-align all the columns of dtype object

    Returns
    -------
    dict : Formatter dictionary

    """
    if cols is None:
       cols = df.columns[df.dtypes == 'object'] 
    return {col: f'{{:<{df[col].str.len().max()}s}}'.format for col in cols}

def main(args):
    table = defaultdict(list)
    aromatic_residues = ["F", "W", "Y"]
    for record in SeqIO.parse(args.fasta, "fasta"):
        table['ID'].append(record.id)
        aromatic_count = sum([record.seq.count(residue) for residue in aromatic_residues])
        table['Aromatics'].append(aromatic_count)
        table['Length'].append(len(record))
        aromatic_freq = aromatic_count / len(record)
        table['Frequency'].append("{:.2%}".format(aromatic_freq))
        gaps = re.split("|".join(aromatic_residues), str(record.seq))
        max_gap = max([len(s) for s in gaps])
        table['Max. Gap'].append(max_gap)
        table['Conductive'].append((aromatic_freq >= args.frequency) and (max_gap <= args.gap))
    df = pd.DataFrame(table)
    print(df.to_string(index=False, justify='center', formatters=make_lalign_formatter(df, cols=['ID'])))

if __name__ == "__main__":
    desc = ("Evaluates input sequences for whether they fit tentative criteria to make an e-pilus. "
            "The criteria used are 1) frequency of aromatic amino acids, and 2) no large gaps "
            "between successive aromatic amino acids.")
    epil = ("Keep in mind that the default threshold values are based on a limited survey from "
            "Walker et al. 2018 (doi:10.1038/ismej.2017.141) and should be treated as conservative "
            "bounds, not exact bounds.")
    parser = argparse.ArgumentParser(description=desc, epilog=epil)
    parser.add_argument("fasta", metavar="FASTA", help="A FASTA file containing one or more sequences.")
    parser.add_argument("-g", "--gap", default=40, type=int, help="Largest allowable gap between aromatic residues. Default: 40")
    parser.add_argument("-f", "--frequency", default=0.09, type=float, help='Minimum allowable frequency of aromatic residues. Default: 0.09 (9%%)')
    args = parser.parse_args()
    if not os.path.isfile(args.fasta):
        sys.exit("The specified file does not exist: {}".format(args.fasta))
    if args.gap < 1:
        sys.exit("The specified gap size must be a positive number.")
    if (args.frequency < 0) or (args.frequency > 1):
        sys.exit("The specified frequency must be between 0 and 1, inclusive.")
    main(args)

