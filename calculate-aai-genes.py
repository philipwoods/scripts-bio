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
    search_results = os.path.join(outdir, name + ".tsv")
    table_path = os.path.join(outdir, name + "-table.tsv")
    # Create the output directory if necessary
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # Create the BLAST database in a temporary directory and query it
    with tempfile.TemporaryDirectory() as tempdir:
        # Construct shell commands for each backend
        if args.backend == 'blastp':
            makedb_args = {
                'seqs': inpath,
                'dbname': name
            }
            makedb_command = ("makeblastdb "
                    "-in {seqs} "
                    "-dbtype prot"
                    "-out {dbname}")
            search_args = {
                'query': inpath,
                'dbname': name,
                'out': search_results,
                'threads': args.threads,
                'maxtargets': args.max_target_seqs
            }
            search_command = ("blastp "
                    "-query {query} "
                    "-db {dbname} "
                    "-out {out} "
                    "-outfmt '6 qseqid sseqid pident' "
                    "-max_hsps 1 "
                    "-num_threads {threads}"
                    "-max_target_seqs {maxtargets}")
        elif args.backend == 'diamond':
            makedb_args = {
                'seqs': inpath,
                'dbname': name,
                'threads': args.threads
            }
            makedb_command = ("diamond makedb "
                    "--in {seqs} "
                    "--db {dbname} "
                    "--threads {threads}")
            search_args = {
                'seqs': inpath,
                'dbname': name,
                'out': search_results,
                'threads': args.threads,
                'maxtargets': args.max_target_seqs
            }
            search_command = ("diamond blastp "
                    "--query {seqs} "
                    "--db {dbname} "
                    "--out {out} "
                    "--threads {threads} "
                    "--sensitive "
                    "--max-hsps 1 "
                    "--outfmt 6 qseqid sseqid pident "
                    "--max-target-seqs {maxtargets}")
        else:
            sys.exit("ERROR: Unsupported backend")
        # Run database and comparison commands
        os.chdir(tempdir)
        subprocess.run(makedb_command.format(**makedb_args), shell=True)
        subprocess.run(search_command.format(**search_args), shell=True)
        os.chdir(outdir)
    tempdf = pd.read_csv(search_results, sep='\t', names=['Query', 'Target', 'Identity'])
    table = tempdf.pivot(index='Query', columns='Target', values='Identity')
    table.to_csv(table_path, sep='\t', float_format="%.3f")

if __name__ == "__main__":
    desc = ("This program accepts a file with multiple homologous protein sequences and computes pairwise amino acid identity between them.")
    epil = ("Note: The different backends have different thresholds below which sequence identity is considered unreliable. "
            "If your input sequences are sufficiently distinct, there may be no reported AAI value for that pairwise comparison. "
            "Missing data for a pairwise comparison may also result if the number of target sequences allowed is smaller than the "
            "number of sequences provided in your input set.")
    parser = argparse.ArgumentParser(description=desc, epilog=epil)
    parser.add_argument("fasta", metavar="FASTA", help="A multi-sequence FASTA file.")
    parser.add_argument("directory", help="The directory to put output files in.")
    parser.add_argument("-b", "--backend", choices=['blastp', 'diamond'], default='diamond', help="The program to use for sequence comparison. Default: %(default)s")
    parser.add_argument("-t", "--threads", default=1, type=int, help="Number of threads for the BLAST search. Default: %(default)s")
    parser.add_argument("--max-target-seqs", default=100, type=int, help="Maximum number of matches allowed per query sequence. Default: %(default)s")
    args = parser.parse_args()
    if not os.path.isfile(args.fasta):
        sys.exit("Specified file does not exist: {}".format(args.fasta))
    if args.backend == 'blastp' and shutil.which("blastp") is None:
        sys.exit("BLAST not found. Please install BLAST or activate an appropriate virtual environment.")
    if args.backend == 'diamond' and shutil.which("diamond") is None:
        sys.exit("DIAMOND not found. Please install DIAMOND or activate an appropriate virtual environment.")
    main(args)

