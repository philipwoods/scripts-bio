#!/usr/bin/env python
import os
import sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
from collections import Counter

def format_columns(table):
    genetic_code = CodonTable.unambiguous_dna_by_id[table].forward_table
    tuples = []
    for codon, aa in genetic_code.items():
        tuples.append((aa, codon))
    return tuples

def split_codons(seq):
    codons = []
    for start_i in range(0, len(seq), 3):
        end_i = start_i + 3
        codons.append(seq[start_i:end_i])
    return codons

def count_codons(codon_set, table):
    genetic_code = CodonTable.unambiguous_dna_by_id[table].forward_table
    c = Counter(codon_set)
    codon_counts = []
    for codon, aa in genetic_code.items():
        codon_counts.append(c[codon])
    return codon_counts

def main(args):
    codon_sets = {}
    for record in SeqIO.parse(args.genes, 'fasta'):
        if record.id in args.rna_genes:
            continue
        codon_sets[record.id] = split_codons(str(record.seq))
        #print(record.id, len(codon_sets[record.id]), codon_sets[record.id][-1])
    codon_counts = {}
    for gene, codons in codon_sets.items():
        codon_counts[gene] = count_codons(codons, args.translation_table)
    data = {
        'index': codon_counts.keys(),
        'index_names': ['gene'],
        'columns': format_columns(args.translation_table),
        'column_names': ['AA', 'codon'],
        'data': codon_counts.values()
    }
    counts_df = pd.DataFrame.from_dict(data, orient='tight')
    counts_df.to_csv(args.output, sep='\t')

if __name__ == "__main__":
    desc = ("")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("genes", help="FASTA file containing coding sequences for each gene to analyze.")
    parser.add_argument("--translation_table", type=int, default=11, help="Translation table to use. Default: %(default)s")
    parser.add_argument("--rna_genes", nargs='+', default=[], help="Gene sequence IDs known to be RNAs, which will be removed from codon analyses.")
    parser.add_argument("--output", default=sys.stdout, help="Location to write output. Default: STDOUT")
    args = parser.parse_args()
    if not os.path.isfile(args.genes):
        sys.exit("The specified file does not exist: {}".format(args.genes))
    main(args)
