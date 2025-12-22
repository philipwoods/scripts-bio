#!/usr/bin/env python3
import os
import sys
import argparse
import tempfile
import subprocess
import pandas as pd
import shutil

def main(args):
    # Get absolute locations for relevant files and directories
    inpath = os.path.abspath(args.fasta)
    outdir = os.path.abspath(args.directory)
    name = os.path.splitext(os.path.basename(inpath))[0]
    blast_results = os.path.join(outdir, name + ".tsv")
    table_path = os.path.join(outdir, name + "-table.tsv")
    # Create the output directory if necessary
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # Create the BLAST database in a temporary directory and query it
    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        subprocess.run("makeblastdb -in {seqs} -dbtype prot -out {dbname}".format(seqs=inpath,dbname=name), shell=True)
        blastp_args = {
            'query': inpath,
            'dbname': name,
            'out': blast_results,
            'threads': args.threads
        }
        subprocess.run("blastp -query {query} -db {dbname} -out {out} -outfmt '6 qseqid sseqid pident' -max_hsps 1 -num_threads {threads}".format(**blastp_args), shell=True)
        os.chdir(outdir)
    tempdf = pd.read_csv(blast_results, sep='\t', names=['Query', 'Target', 'Identity'])
    table = tempdf.pivot(index='Query', columns='Target', values='Identity')
    table.to_csv(table_path, sep='\t', float_format="%.3f")

if __name__ == "__main__":
    desc = ("")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("fasta", metavar="FASTA", help="A multi-sequence FASTA file.")
    parser.add_argument("directory", help="The directory to put output files in.")
    parser.add_argument("-t", "--threads", default=1, type=int, help="Number of threads for the BLAST search. Default: 1")
    args = parser.parse_args()
    if not os.path.isfile(args.fasta):
        sys.exit("Specified file does not exist: {}".format(args.fasta))
    if shutil.which("blastp") is None:
        sys.exit("BLAST not found. Please install BLAST or activate an appropriate virtual environment.")
    main(args)

