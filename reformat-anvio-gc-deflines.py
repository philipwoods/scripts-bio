#!/usr/bin/env python3
import sys
import os
import argparse

# For the record, the defline format from anvi-get-sequences-for-gene-clusters
# is the following when you ask for unconcatenated alignments:
# >TEMPID|gene_cluster:GC|genome_name:GENOME_NAME|gene_callers_id:GENE_ID
# This is not super useful sometimes, so we want to convert to:
# >GENOME_NAME___GENE_ID|gene_cluster:GC

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
                fields = line.split("|")
                genome = fields[2].split(":")[1].strip()
                geneid = fields[3].split(":")[1].strip()
                cluster = fields[1].strip()
                output.write(">{name}___{gene}|{gc}\n".format(name=genome,gene=geneid,gc=cluster))
    if args.inplace:
        os.replace(tmpfile, args.fasta)
        output.close()

if __name__ == "__main__":
    desc = ("This script reformats the deflines in a FASTA from anvi-get-sequences-for-gene-clusters "
            "(without the --concatenate flag) to be more useful than the default.")
    epil = ("The default defline format from anvi-get-sequences-for-gene-clusters (without the "
            "--concatenate flag) is: '>######|gene_cluster:GC|genome_name:GENOME_NAME|gene_callers_id:GENE_ID'. "
            "This script reformats these to: '>GENOME_NAME___GENE_ID|gene_cluster:GC'")
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

