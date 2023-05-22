#!/usr/bin/env python3

import os
import sys
import argparse
import textwrap
from Bio import SeqIO

format_exts = {
    'clustal': set(['.clustal', '.aln']),
    'embl': set(['.embl']),
    'fasta': set(['.fa', '.faa', '.fna', '.fasta']),
    'fasta-2line': set([]),
    'genbank': set(['.gb', '.genbank']),
    'imgt': set(['.imgt']),
    'maf': set(['.maf']),
    'mauve': set(['.alignment', '.xmfa']),
    'nexus': set(['.nexus']),
    'phylip': set(['.phy', '.ph', '.phylip']),
    'phylip-sequential': set([]),
    'phylip-relaxed': set([]),
    'pir': set(['.pir', '.nbrf']),
    'seqxml': set(['.xml']),
    'stockholm': set(['.sto', '.stk']),
    'tab': set(['.tab', '.tsv']),
    'xdna': set(['.xdna'])
}

def infer_format(filename):
    extension = os.path.splitext(filename)[1]
    for fmt, exts in format_exts.items():
        if extension in exts:
            return fmt
    return None

def main(args):
    if args.infmt is None:
        fmt = infer_format(args.input)
        if fmt is None:
            sys.exit("Input file format not recognized: {}".format(args.input))
        args.infmt = fmt
    if args.outfmt is None:
        if args.output == sys.stdout:
            sys.exit("You must specify either an output file destination or an output format.")
        fmt = infer_format(args.output.name)
        if fmt is None:
            sys.exit("Output file format not recognized: {}".format(args.output.name))
        args.outfmt = fmt
    SeqIO.convert(args.input, args.infmt, args.output, args.outfmt)

if __name__ == "__main__":
    desc =  textwrap.dedent("""\
        This script converts between different sequence file formats. It will infer the
        input and desired output file formats from the extensions of the provided file names.
        Users can also specify file formats, in which case the file extensions are ignored.
        """)
    epil = textwrap.dedent("""\
        Supported format strings and inferred file extensions
        -----------------------------------------------------------
            clustal             {clustal}
            embl                {embl}
            fasta               {fasta}
            fasta-2line         (Note: exactly 2 lines per record, no wrapping)
            genbank             {genbank}
            imgt                {imgt}
            maf                 {maf}
            mauve               {mauve}
            nexus               {nexus}
            phylip              {phylip}
            phylip-sequential
            phylip-relaxed      (Note: interleaved, but allows longer sequence names)
            pir                 {pir}
            seqxml              {seqxml}
            stockholm           {stockholm}
            tab                 {tab}
            xdna                {xdna}
        """.format(**format_exts))
    parser = argparse.ArgumentParser(description=desc, epilog=epil, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--infmt', metavar="format-string", choices=format_exts.keys(), help="")
    parser.add_argument('--outfmt', metavar="format-string", choices=format_exts.keys(), help="")
    parser.add_argument('input', help="A file containing one or more sequences or sequence alignments.")
    parser.add_argument('output', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="A reformatted file containing the same sequence(s).")
    args = parser.parse_args()
    main(args)

