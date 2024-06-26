#!/usr/bin/env python3
# Import required modules
import os
import sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import pandas as pd
from collections import defaultdict

def main():
    help_string = """
    Synopsis:
        python parse-alignment-column-frequencies.py <output-dir> <alignment> <groups-file>
    
    Description:
        This is an auxiliary utility for the script anvi-script-alignment-enrichment.sh
        and is not generally intended for independent use. For example, it has no input
        validation and relies on its parent script always passing well-formed input. For
        that script to work, both scripts must be in the same directory.
    """
    version_string = """
    Last updated 19 June 2024
    """

    if "-h" in sys.argv or "--help" in sys.argv:
        print(help_string)
        sys.exit(1)

    #Handle arguments
    output_dir = sys.argv[1]
    alignment_file= sys.argv[2]
    groups_file = sys.argv[3]

    # Declare amino acid reference dictionary for output
    amino_acids = {
        'key': ['Alanine', 'Cysteine', 'Aspartate', 'Glutamate', 'Phenylalanine',
                'Glycine', 'Histidine', 'Isoleucine', 'Lysine', 'Leucine',
                'Methionine', 'Asparagine', 'Proline', 'Glutamine', 'Arginine',
                'Serine', 'Threonine', 'Valine', 'Tryptophan', 'Tyrosine', 'Gap'],
        'accession': ['A', 'C', 'D', 'E', 'F',
                      'G', 'H', 'I', 'K', 'L',
                      'M', 'N', 'P', 'Q', 'R',
                      'S', 'T', 'V', 'W', 'Y', '-']
    }

    # Read in the alignment and sort by the ID field
    with open(alignment_file, 'r') as f:
        alignment = AlignIO.read(f, 'fasta')
    # Read in groups file, columns 'samples' and 'group'
    groups_db = pd.read_csv(groups_file, sep='\t', index_col='samples')
    # Get a dictionary of the form {'G1': [seqid0, seqid2,...], 'G2': [seqid1, seqid3,...]}
    groups = groups_db.groupby('group').groups

    # create subset alignments of each group
    group_seqs = defaultdict(list)
    group_alns = {}
    for seqrecord in alignment:
        for group in groups.keys():
            if seqrecord.id in groups[group]:
                group_seqs[group].append(seqrecord)
    print(group_seqs)
    for group in groups.keys():
        group_alns[group] = MultipleSeqAlignment(group_seqs[group])
        print(group_alns[group])

    out_df = pd.DataFrame(data=amino_acids)
    out_df['function'] = "NA"

if __name__ == "__main__":
    main()
