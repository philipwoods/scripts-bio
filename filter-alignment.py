#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
from Bio import AlignIO

# Amino acid conservation properties extracted from jalview/schemes/ResidueProperties.java
jalview_props = ["hydrophobic", "polar", "small", "positive", "negative", "charged", "aromatic", "aliphatic", "tiny", "proline"]
aa_props = {
    "A": set(["hydrophobic", "not polar", "small", "not positive", "not negative", "not charged", "not aromatic", "not aliphatic", "tiny", "not proline"]),
    "C": set(["hydrophobic", "not polar", "small", "not positive", "not negative", "not charged", "not aromatic", "not aliphatic", "not tiny", "not proline"]),
    "D": set(["not hydrophobic", "polar", "small", "not positive", "negative", "charged", "not aromatic", "not aliphatic", "not tiny", "not proline"]),
    "E": set(["not hydrophobic", "polar", "not small", "not positive", "negative", "charged", "not aromatic", "not aliphatic", "not tiny", "not proline"]),
    "F": set(["hydrophobic", "not polar", "not small", "not positive", "not negative", "not charged", "aromatic", "not aliphatic", "not tiny", "not proline"]),
    "G": set(["hydrophobic", "not polar", "small", "not positive", "not negative", "not charged", "not aromatic", "not aliphatic", "tiny", "not proline"]),
    "H": set(["hydrophobic", "polar", "not small", "positive", "not negative", "charged", "aromatic", "not aliphatic", "not tiny", "not proline"]),
    "I": set(["hydrophobic", "not polar", "not small", "not positive", "not negative", "not charged", "not aromatic", "aliphatic", "not tiny", "not proline"]),
    "K": set(["hydrophobic", "polar", "not small", "positive", "not negative", "charged", "not aromatic", "not aliphatic", "not tiny", "not proline"]),
    "L": set(["hydrophobic", "not polar", "not small", "not positive", "not negative", "not charged", "not aromatic", "aliphatic", "not tiny", "not proline"]),
    "M": set(["hydrophobic", "not polar", "not small", "not positive", "not negative", "not charged", "not aromatic", "not aliphatic", "not tiny", "not proline"]),
    "N": set(["not hydrophobic", "polar", "small", "not positive", "not negative", "not charged", "not aromatic", "not aliphatic", "not tiny", "not proline"]),
    "P": set(["not hydrophobic", "not polar", "small", "not positive", "not negative", "not charged", "not aromatic", "not aliphatic", "not tiny", "proline"]),
    "Q": set(["not hydrophobic", "polar", "not small", "not positive", "not negative", "not charged", "not aromatic", "not aliphatic", "not tiny", "not proline"]),
    "R": set(["not hydrophobic", "polar", "not small", "positive", "not negative", "charged", "not aromatic", "not aliphatic", "not tiny", "not proline"]),
    "S": set(["not hydrophobic", "polar", "small", "not positive", "not negative", "not charged", "not aromatic", "not aliphatic", "tiny", "not proline"]),
    "T": set(["not hydrophobic", "polar", "small", "not positive", "not negative", "not charged", "not aromatic", "not aliphatic", "not tiny", "not proline"]),
    "V": set(["hydrophobic", "not polar", "small", "not positive", "not negative", "not charged", "not aromatic", "aliphatic", "not tiny", "not proline"]),
    "W": set(["hydrophobic", "polar", "not small", "not positive", "not negative", "not charged", "aromatic", "not aliphatic", "not tiny", "not proline"]),
    "Y": set(["hydrophobic", "polar", "not small", "not positive", "not negative", "not charged", "aromatic", "not aliphatic", "not tiny", "not proline"])
    }

def main(args):
    # Read in alignment file
    alignment = AlignIO.read(args.fasta, "fasta")
    ## Late argument validation
    # If no outgroup is given, the outgroup should be empty.
    if args.outgroup is None:
        args.outgroup = []
    out_set = set(args.outgroup)
    # If no ingroup is given, the ingroup should be everything not in the outgroup.
    if args.ingroup is None:
        in_set = set(range(len(alignment)))
        in_set = in_set.difference(out_set)
        args.ingroup = sorted(in_set)
    in_set = set(args.ingroup)
    if not in_set.isdisjoint(out_set):
        sys.exit("Sequences cannot be in both the ingroup and the outgroup.")
    # Handle list display argument
    if args.list:
        for i, record in enumerate(alignment):
            print("{index}\t{header}".format(index=i, header=record.id))
        sys.exit()
    ## Compute identity at each site in the alignment
    conserved_sites = []
    for i in range(alignment.get_alignment_length()):
        # Extract one column as a string
        column = alignment[:,i]
        # Calculate the frequency of gaps in the column
        # and disregard it if there are too many
        gap_freq = column.count('-') / len(column)
        if gap_freq > args.gaps:
            continue
        # Identify the unique residues in the column
        residues = set(column)
        if '-' in residues:
            residues.remove('-')
        # Get the frequency of each amino acid in the column
        freqs = {}
        for aa in sorted(residues):
            freqs[aa] = column.count(aa) / len(column)
        # Disregard the column if max frequency is less than identity cutoff
        if max(freqs.values()) < args.identity:
            continue
        # If there is a conserved residue, identify it
        # Break statement ensures only one residue per site in case of a tie
        for aa, freq in freqs.items():
            if freq == max(freqs.values()):
                conserved_sites.append({'position': i+1, 'residue': aa, 'frequency': freq})
                break
    # Print conserved sites
    if not args.print_alignment:
        print('Position\tConserved residue\tFrequency')
        for site in conserved_sites:
            print("{position}\t{aa}\t{freq:.2f}".format(position=site['position'], aa=site['residue'], freq=site['frequency']))
    else:
        trimmed = alignment[:,0:1] # Create initial placeholder slice to add to
        positions = []
        consensus = []
        frequencies = []
        # Build trimmed alignment of only conserved sites
        for site in conserved_sites:
            trimmed = trimmed + alignment[:,site['position']-1:site['position']]
            positions.append(site['position'])
            consensus.append(site['residue'])
            frequencies.append(site['frequency'])
        trimmed = trimmed[:,1:] # Remove initial placeholder column
        # Handle annotations and display
        annotations = pd.DataFrame({'Position': positions, 'Consensus': consensus, 'Frequency': frequencies})
        print(format(trimmed, 'clustal'))
        print("Consensus sequence:\t{}\n".format("".join(consensus)))
        print(annotations.T.to_string(header=False, float_format='{:.2f}'.format))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', metavar='FASTA', help="FASTA file containing the alignment of interest")
    parser.add_argument('--list', '-l', action='store_true', help="Display indices of sequence headers in the FASTA")
    parser.add_argument('--ingroup', type=int, nargs='+', help="Specify indices of sequences to check for site conservation. Default: all")
    parser.add_argument('--outgroup', type=int, nargs='+', help="Specify indices of sequences which should not be conserved with the ingroup. Default: none")
    parser.add_argument('--identity', '-i', type=float, default=1, help="Set the minimum identity threshold for conservation. Default: %(default)s")
    parser.add_argument('--gaps', '-g', type=float, default=0.25, help="Set the maximum proportion of gaps allowed in a site for conservation analysis. Default: %(default)s")
    parser.add_argument('--print-alignment', '-p', action='store_true', help="Print a trimmed alignment of the conserved sites instead of a summary table.")
    args = parser.parse_args()
    # Argument validation
    if not os.path.isfile(args.fasta):
        sys.exit("Specified file does not exist: {}".format(args.fasta))
    if (args.identity < 0) or (args.identity > 1):
        sys.exit("Identity threshold must be a proportion between 0 and 1, inclusive.")
    if (args.gaps < 0) or (args.gaps > 1):
        sys.exit("Gaps threshold must be a proportion between 0 and 1, inclusive.")
    main(args)

