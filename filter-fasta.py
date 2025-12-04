#!/usr/bin/env python3
import sys
import os.path
import argparse
import numpy as np
from Bio import SeqIO

class Histogram(object):
    """
    Ascii histogram
    Taken from https://pyinsci.blogspot.com/2009/10/ascii-histograms.html
    Licenced under GPL
    """
    def __init__(self, data, bins=10):
        """
        Class constructor

        :Parameters:
            - `data`: array like object
        """
        self.data = data
        self.bins = bins
        self.h = np.histogram(self.data, bins=self.bins)

    def horizontal(self, height=4, character ='|'):
        """Returns a multiline string containing a
        a horizontal histogram representation of self.data
        :Parameters:
            - `height`: Height of the histogram in characters
            - `character`: Character to use
        >>> d = normal(size=1000)
        >>> h = Histogram(d,bins=25)
        >>> print h.horizontal(5,'|')
        106            |||
                      |||||
                      |||||||
                    ||||||||||
                   |||||||||||||
        -3.42                         3.09
        """
        his = """"""
        bars = self.h[0]*height/max(self.h[0])
        for l in reversed(range(1,height+1)):
            line = ""
            if l == height:
                line = '%s '%max(self.h[0]) #histogram top count
            else:
                line = ' '*(len(str(max(self.h[0])))+1) #add leading spaces
            for c in bars:
                if c >= np.ceil(l):
                    line += character
                else:
                    line += ' '
            line +='\n'
            his += line
        his += '%.2f'%self.h[1][0] + ' '*(self.bins) +'%.2f'%self.h[1][-1] + '\n'
        return his

    def vertical(self,height=20, character ='|'):
        """
        Returns a Multi-line string containing a
        a vertical histogram representation of self.data
        :Parameters:
            - `height`: Height of the histogram in characters
            - `character`: Character to use
        >>> d = normal(size=1000)
        >>> Histogram(d,bins=10)
        >>> print h.vertical(15,'*')
                              236
        -3.42:
        -2.78:
        -2.14: ***
        -1.51: *********
        -0.87: *************
        -0.23: ***************
        0.41 : ***********
        1.04 : ********
        1.68 : *
        2.32 :
        """
        his = """"""
        xl = ['%.2f'%n for n in self.h[1]]
        lxl = [len(l) for l in xl]
        bars = self.h[0]*height//max(self.h[0])
        his += ' '*(max(bars)+2+max(lxl))+'%s\n'%max(self.h[0])
        for i,c in enumerate(bars):
            line = xl[i] +' '*(max(lxl)-lxl[i])+': '+ character*c+'\n'
            his += line
        return his

def within_length_bounds(record, args):
    return (len(record) > args.min) and (len(record) < args.max)

def main(args):
    if args.list:
        for record in SeqIO.parse(args.fasta, "fasta"):
            print(record.id, len(record), sep='\t')
        sys.exit()
    if args.summarize:
        lengths = []
        for record in SeqIO.parse(args.fasta, "fasta"):
            lengths.append(len(record))
        lenrange = int(np.max(lengths) - np.min(lengths))
        h = Histogram(lengths, bins=30)
        print(h.vertical(120))
        print("")
        print("Min:    {: 4d}".format(np.min(lengths)))
        print("Max:    {: 4d}".format(np.max(lengths)))
        print("Mean:   {: 1g}".format(np.mean(lengths)))
        print("Stdev:  {: 1g}".format(np.std(lengths)))
        print("Median: {: 1g}".format(np.median(lengths)))
        print("P10:    {: 1g}".format(np.percentile(lengths, 10)))
        print("P25:    {: 1g}".format(np.percentile(lengths, 25)))
        print("P75:    {: 1g}".format(np.percentile(lengths, 75)))
        print("P90:    {: 1g}".format(np.percentile(lengths, 90)))
        sys.exit()
    # Add entries from the exclusion file to the excluded ID list
    if args.exclude_file is not None:
        with open(args.exclude_file) as f:
            file_ids = [line.strip() for line in f.readlines()]
        if args.exclude_ids is None:
            args.exclude_ids = file_ids
        else:
            args.exclude_ids.extend(file_ids)
    # Add entries from the inclusion file to the included ID list
    if args.include_file is not None:
        with open(args.include_file) as f:
            file_ids = [line.strip() for line in f.readlines()]
        if args.include_ids is None:
            args.include_ids = file_ids
        else:
            args.include_ids.extend(file_ids)
    # Create list of IDs to output
    # Note that argument checking excludes the case where both args.include_ids and args.exclude_ids have a value
    records = SeqIO.parse(args.fasta, "fasta")
    if args.include_ids is None and args.exclude_ids is None:
        output_ids = [record.id for record in records]
    elif args.include_ids is not None:
        output_ids = [record.id for record in records if record.id in args.include_ids]
    elif args.exclude_ids is not None:
        output_ids = [record.id for record in records if record.id not in args.exclude_ids]
    output_records = []
    # We need to create a new iterator after running through it creating output_ids
    for record in SeqIO.parse(args.fasta, "fasta"):
        # The record ID is anything before a space in the defline
        if (record.id in output_ids) and within_length_bounds(record, args):
            output_records.append(record)
    SeqIO.write(output_records, args.out, "fasta")

if __name__ == "__main__":
    desc = ("Takes in a multi-sequence FASTA file and searches it for the provided sequence IDs "
            "and/or by the provided minimum and maximum sequence lengths.")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("fasta", metavar="FASTA", help="A multi-sequence FASTA file.")
    parser.add_argument("--list", "-l", action='store_true', help="List the available sequence IDs from the input.")
    parser.add_argument("--summarize", action='store_true', help="Summarize length statistics with a histogram.")
    parser.add_argument("--out", default=sys.stdout, help="A path to output the results. Default: stdout")
    parser.add_argument("--min", type=int, default=0, help="Minimum acceptable sequence length to output.")
    parser.add_argument("--max", type=int, default=100000000000000000, help="Maximum acceptable sequence length to output.")
    group_include = parser.add_argument_group(title="Inclusion")
    group_include.add_argument("--include-ids", nargs='+', help="A list of sequence IDs to output. Default: output everything")
    group_include.add_argument("--include-file", help="A file containing a list of sequence IDs, one per line, to output.")
    group_exclude = parser.add_argument_group(title="Exclusion")
    group_exclude.add_argument("--exclude-ids", nargs='+', help="A list of sequence IDs to exclude from output. Default: output everything")
    group_exclude.add_argument("--exclude-file", help="A file containing a list of sequence IDs, one per line, to exclude from output.")
    args = parser.parse_args()
    if not os.path.isfile(args.fasta):
        sys.exit("The specified file does not exist: {}".format(args.fasta))
    if args.include_file is not None and not os.path.isfile(args.include_file):
        sys.exit("The specified file does not exist: {}".format(args.include_file))
    if args.exclude_file is not None and not os.path.isfile(args.exclude_file):
        sys.exit("The specified file does not exist: {}".format(args.exclude_file))
    if (args.include_ids or args.include_file) and (args.exclude_ids or args.exclude_file):
        sys.exit("You must specify either inclusion or exclusion criteria, but not both.")
    main(args)

