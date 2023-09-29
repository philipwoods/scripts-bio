#!/usr/bin/env python3
import sys
import os
import argparse
import re

# For the record, the defline format from anvi-get-sequences-for-hmm-hits
# is the following when you ask for unconcatenated sequences:
# >GENENAME___HMMSOURCE___LONGHASH bin_id:GENOMENAME|source:HMMSOURCE|e_value:XXXX|contig:HASH_CONTIGNAME|gene_callers_id:HASH_GENEID|start:START|stop:STOP|length:YYYY
# This is very long and the space in the name cuts of relevant info, so we want to convert to:
# >GENOMENAME___HMMSOURCE___GENENAME|gene_id:HASH_GENEID

def main(args):
    if args.inplace:
        tmpfile = args.fasta + ".tmp"
        output = open(tmpfile, 'w')
    else:
        output = args.output
    with open(args.fasta, 'r') as f:
        for line in f:
            if line[0] != ">":
                output.write(line)
            else:
                line = line[1:]
                fields = re.split(r'___| |\|',line) # Split on '|' or ' ' or '___'
                genename= fields[0]
                hmm = fields[1]
                genome = fields[3].split(":")[1].strip()
                geneid = fields[7].split(":")[1].strip()
                output.write(">{name}___{hmm}___{gene}|gene_id:{geneid}\n".format(name=genome,hmm=hmm,gene=genename,geneid=geneid))
    if args.inplace:
        os.replace(tmpfile, args.fasta)
        output.close()

if __name__ == "__main__":
    desc = ("This script reformats the deflines in a FASTA from anvi-get-sequences-for-hmm-hits "
            "(without the --concatenate flag) to be more useful than the default for phylogeny.")
    epil = "This script reformats deflines to: '>GENOMENAME___HMMCSOURCE___GENE|gene_id:GENEID'"
    parser = argparse.ArgumentParser(description=desc, epilog=epil)
    parser.add_argument('fasta', help="FASTA file from anvi-get-sequences-from-gene-clusters")
    parser.add_argument('-i', '--inplace', action='store_true',
                        help="Modify the input file in-place rather than writing to output.")
    parser.add_argument('-o', '--output', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Write converted FASTA to this file. Default: STDOUT")
    args = parser.parse_args()
    if not os.path.isfile(args.fasta):
        sys.exit("Specified file does not exist: {}".format(args.fasta))
    main(args)

