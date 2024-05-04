#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

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

# From Susko & Roger (2007) doi:10.1093/molbev/msm144
recoding_bins = {
        "SR4": {"A": "1", "C": "2", "D": "3", "E": "3", "F": "4", "G": "1", "H": "2", "I": "4", "K": "3", "L": "4", "M": "4", "N": "1", "P": "1", "Q": "3", "R": "3", "S": "1", "T": "1", "V": "4", "W": "2", "Y": "2"},
        "SR5": {"A": "1", "C": "2", "D": "3", "E": "3", "F": "2", "G": "1", "H": "4", "I": "5", "K": "4", "L": "5", "M": "5", "N": "3", "P": "1", "Q": "4", "R": "4", "S": "1", "T": "1", "V": "5", "W": "2", "Y": "2"},
        "SR6": {"A": "1", "C": "2", "D": "3", "E": "3", "F": "4", "G": "3", "H": "4", "I": "5", "K": "6", "L": "5", "M": "5", "N": "3", "P": "1", "Q": "6", "R": "6", "S": "1", "T": "1", "V": "5", "W": "2", "Y": "4"},
        "SR7": {"A": "1", "C": "2", "D": "3", "E": "3", "F": "4", "G": "1", "H": "5", "I": "6", "K": "7", "L": "6", "M": "6", "N": "3", "P": "5", "Q": "7", "R": "7", "S": "1", "T": "1", "V": "6", "W": "2", "Y": "4"},
        "SR8": {"A": "1", "C": "2", "D": "3", "E": "3", "F": "4", "G": "2", "H": "5", "I": "6", "K": "7", "L": "6", "M": "8", "N": "3", "P": "5", "Q": "7", "R": "7", "S": "1", "T": "1", "V": "6", "W": "8", "Y": "4"},
        "SR9": {"A": "1", "C": "2", "D": "3", "E": "3", "F": "4", "G": "5", "H": "6", "I": "7", "K": "8", "L": "7", "M": "9", "N": "5", "P": "9", "Q": "6", "R": "8", "S": "1", "T": "1", "V": "7", "W": "2", "Y": "4"},
        "SR10": {"A": "1", "C": "2", "D": "3", "E": "3", "F", "4", "G", "5", "H", "6", "I": "7", "K": "8", "L": "9", "M": "9", "N": "5", "P": "P", "Q": "6", "R": "8", "S": "1", "T": "1", "V": "7", "W": "2", "Y": "4"},
        "SR11": {"A": "1", "C": "C", "D": "2", "E": "2", "F": "3", "G": "4", "H": "5", "I": "6", "K": "7", "L": "8", "M": "8", "N": "4", "P": "P", "Q": "5", "R": "7", "S": "1", "T": "1", "V": "6", "W": "W", "Y": "3"},
        "SR13": {"A": "1", "C": "C", "D": "2", "E": "2", "F": "3", "G": "G", "H": "H", "I": "4", "K": "5", "L": "6", "M": "6", "N": "N", "P": "P", "Q": "Q", "R": "5", "S": "1", "T": "1", "V": "4", "W": "W", "Y": "3"},
        "SR15": {"A": "1", "C": "C", "D": "2", "E": "2", "F": "F", "G": "G", "H": "H", "I": "3", "K": "4", "L": "L", "M": "M", "N": "N", "P": "P", "Q": "Q", "R": "4", "S": "1", "T": "1", "V": "3", "W": "W", "Y": "Y"}
    }

def residue_freqs(column):
    """
    Arguments
    ----------------
    column : str
        A string containing the residues in one column of an alignment

    Returns
    ----------------
    A dict with residues as keys and frequencies as values
    """
    # Identify the unique residues in the column
    residues = set(column)
    # Get the frequency of each amino acid in the column
    freqs = {}
    for aa in sorted(residues):
        freqs[aa] = column.count(aa) / len(column)
    return freqs

def conserved_properties(column, max_gaps):
    """
    Arguments
    ----------------
    column : str
        A string containing the residues in one column of an alignment
    max_gaps : float
        The maximum proportion of gaps allowed in the column

    Returns
    ----------------
    Tuple containing:
        1) A float indicating the proportion of chemical properties conserved in the column
        2) A one-character string indicating the degree of chemical property conservation
    """
    if len(column) < 1:
        return

    # Calculate the frequency of gaps in the column
    # and disregard it if there are too many
    freqs = residue_freqs(column) 
    if ('-' in freqs) and (freqs['-'] > max_gaps):
        return (0, '-')

    # Identify the unique residues in the column
    residues = set(column)
    # Remove gaps or undetermined residues
    residues.discard('-')
    residues.discard('X')
    # If the column is all gaps or all undetermined,
    # residues will now be empty. Handle this as above.
    if len(residues) == 0:
        return (0, '-')

    # Identify shared chemical properties
    conserved = None
    for residue in residues:
        if conserved is None:
            conserved = aa_props[residue]
        else:
            conserved = conserved.intersection(aa_props[residue])
    
    # Return a single-character string
    # The maximum value is 10, encoded as '*'
    propstr = str(len(conserved))
    if len(conserved) > 9:
        propstr = '*'
    return (len(conserved)/len(jalview_props), propstr)

def find_conserved_sites(groups, args):
    """
    Arguments
    ----------------
    groups : dict of Biopython MSA objects containing
        'ingroup':  A MSA object of the taxa to be compared for conservation
        'outgroup': A MSA object of the taxa which cannot share the conserved residue
        'others':   A MSA object of taxa in neither the ingroup or outgroup
    args : ArgParse object
        An object containing the command-line parameters provided by the user

    Returns
    ----------------
    cons_sites : list of conserved sites where each site is encoded as a dict containing 
        'position':   its position in the original alignment,
        'residue':    the consensus residue, and
        'frequency':  the residue's frequency in the ingroup.
        'properties': the degree of conserved properties at this site
    """
    ingroup = groups['ingroup']
    outgroup = groups['outgroup']
    conserved_sites = []
    for i in range(ingroup.get_alignment_length()):
        # Extract one column as a string
        column = str(ingroup[:,i])
        # Calculate the frequency of each residue in the column
        freqs = residue_freqs(column)
        # Disregard the column if there are too many gaps
        if ('-' in freqs) and (freqs['-'] > args.gaps):
            continue
        # Disregard the column if max frequency is less than identity cutoff
        if max(freqs.values()) < args.identity:
            continue
        # Disregard the column if property conservation is too low
        if not args.recode:
            prop_cons = conserved_properties(column, args.gaps)
            if prop_cons[0] < args.properties:
                continue
        # If there is a conserved residue, record and annotate it
        # Break statement ensures only one residue per site in case of a tie
        for aa, freq in freqs.items():
            if (freq == max(freqs.values())) and (aa not in outgroup[:,i]):
                properties = "-"
                if args.recode is None:
                    properties = conserved_properties(column, args.gaps)[1]
                site = {
                        'position': i+1,
                        'residue': aa,
                        'frequency': freq,
                        'properties': properties
                    }
                conserved_sites.append(site)
                break
    return conserved_sites

## should properties conservation replace AAI when used or should both be always considered?

def trim_clustal(alignment, header, footer=None):
    """
    Arguments
    ----------------
    alignment : str
        The string output of format(MSA, 'clustal')
    header : str
        A string which will act as a header for the alignment
    footer : str
        An optional string to act as a footer for the alignment

    Returns
    ----------------
    An adjusted clustal formatted alignment string
    """
    lines = alignment.split('\n')
    lines = lines[1:-2]
    lines[0] = header
    if footer is not None:
        lines.extend([footer, ''])
    return "\n".join(lines)

def format_conserved_alignment(groups, cons_sites, args):
    """
    Arguments
    ----------------
    groups : dict of Biopython MSA objects containing
        'ingroup':  A MSA object of the taxa to be compared for conservation
        'outgroup': A MSA object of the taxa which cannot share the conserved residue
        'others':   A MSA object of taxa in neither the ingroup or outgroup
    cons_sites : list of conserved sites where each site is encoded as a dict containing 
        'position':   its position in the original alignment,
        'residue':    the consensus residue, and
        'frequency':  the residue's frequency in the ingroup.
        'properties': the degree of conserved properties at this site
    args : ArgParse object
        An object containing the command-line parameters provided by the user

    Returns
    ----------------
    None (prints output)
    """
    ingroup = groups['ingroup']
    outgroup = groups['outgroup']
    others = groups['others']
    # Create initial placeholder slices for each group
    trimmed_in = ingroup[:,0:1]
    trimmed_out = outgroup[:,0:1]
    trimmed_other = others[:,0:1]
    positions = []
    consensus = []
    frequencies = []
    properties = {'ingroup': [], 'outgroup': [], 'others': []}
    # Build trimmed alignment of only conserved sites
    for site in cons_sites:
        new_in = ingroup[:,site['position']-1:site['position']]
        new_out = outgroup[:,site['position']-1:site['position']]
        new_other = others[:,site['position']-1:site['position']]
        trimmed_in = trimmed_in + new_in
        trimmed_out = trimmed_out + new_out
        trimmed_other = trimmed_other + new_other
        positions.append(site['position'])
        consensus.append(site['residue'])
        frequencies.append(site['frequency'])
        properties['ingroup'].append(site['properties'])
        if args.recode is None:
            if len(new_out) > 1:
                properties['outgroup'].append(conserved_properties(str(new_out[:,0]),args.gaps)[1])
            if len(new_other) > 1:
                properties['others'].append(conserved_properties(str(new_other[:,0]),args.gaps)[1])
        else:
            if len(new_out) > 1:
                properties['outgroup'].append('-')
            if len(new_other) > 1:
                properties['others'].append('-')
    # Remove initial placeholder columns
    trimmed_in = trimmed_in[:,1:]
    trimmed_out = trimmed_out[:,1:]
    trimmed_other = trimmed_other[:,1:]
    # Handle annotations and display
    if args.fasta_out:
        for record in trimmed_in:
            record.id = "Ingroup_" + record.id
        for record in trimmed_out:
            record.id = "Outgroup_" + record.id
        for record in trimmed_other:
            record.id = "Other_" + record.id
        print(format(trimmed_in, 'fasta'))
        print(">Ingroup_Consensus")
        print("".join(consensus))
        print(format(trimmed_out, 'fasta'))
        print(format(trimmed_other, 'fasta'))
    else:
        annotations = pd.DataFrame({'Position': positions, 'Consensus': consensus, 'Frequency': frequencies, 'Properties': properties['ingroup']})
        cons_footer = "Consensus sequence:                 {}".format("".join(consensus))
        prop_footer = "Conserved properties:               {}".format("".join(properties['ingroup']))
        footer = "\n".join([cons_footer, prop_footer])
        print(trim_clustal(format(trimmed_in, 'clustal'), "Ingroup conserved sites", footer))
        if args.outgroup:
            footer = "Conserved properties:               {}".format("".join(properties['outgroup']))
            print(trim_clustal(format(trimmed_out, 'clustal'), "Outgroup sequences", footer))
        if others:
            footer = "Conserved properties:               {}".format("".join(properties['others']))
            print(trim_clustal(format(trimmed_other, 'clustal'), "Other sequences", footer))
        print(annotations.T.to_string(header=False, float_format='{:.2f}'.format))
        print("")
        print("Settings")
        print("---------------------------")
        print("Residue identity:   >={}".format(args.identity))
        print("Prop. conservation: >={}".format(args.properties))
        print("Max. allowed gaps:    {}".format(args.gaps))

def main(args):
    # Read in alignment file
    alignment = AlignIO.read(args.fasta, "fasta")
    # Handle list display argument
    if args.list:
        for i, record in enumerate(alignment):
            print("{index}\t{header}".format(index=i, header=record.id))
        sys.exit()
    # Recode if necessary
    if args.recode is not None:
        recode_dict = recoding_bins[args.recode]
        for record in alignment:
            newseq = record.seq
            for k,v in recode_dict.items():
                newseq = newseq.replace(k,v)
            record.seq = newseq
        recoding_file = os.path.splitext(args.fasta)[0] + "-recoded-{}.fa".format(args.recode)
        AlignIO.write(alignment, recoding_file, "fasta")

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
    # The ingroup and outgroup cannot overlap.
    if not in_set.isdisjoint(out_set):
        sys.exit("Sequences cannot be in both the ingroup and the outgroup.")
    # Create alignment objects for ingroup and outgroup
    others_list = sorted(set(range(len(alignment))).difference(out_set).difference(in_set))
    ingroup = MultipleSeqAlignment([alignment[i] for i in args.ingroup])
    outgroup = MultipleSeqAlignment([alignment[i] for i in args.outgroup]) 
    others = MultipleSeqAlignment([alignment[i] for i in others_list])
    groups = {'ingroup': ingroup, 'outgroup': outgroup, 'others': others}

    ## Compute identity at each site in the alignment
    # Output has format [{'position': X, 'residue': X, 'frequency': F}, {...}, ...]
    conserved_sites = find_conserved_sites(groups, args)

    ## Print conserved sites
    # Handle case with no matching conserved sites
    if len(conserved_sites) < 1:
        sys.exit("No matching conserved sites were identified.")
    # Print formatted output
    if args.table:
        print('Position\tConserved residue\tFrequency\tProperties')
        for site in conserved_sites:
            print("{pos}\t{aa}\t{freq:.3f}\t{prop}".format(
                    pos=site['position'],
                    aa=site['residue'],
                    freq=site['frequency'],
                    prop=site['properties']
                ))
    else:
        format_conserved_alignment(groups, conserved_sites, args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', metavar='FASTA', help="FASTA file containing the alignment of interest")
    parser.add_argument('--list', '-l', action='store_true', help="Display indices of sequence headers in the FASTA")
    parser.add_argument('--ingroup', type=int, nargs='+', help="Specify indices of sequences to check for site conservation. Default: all")
    parser.add_argument('--outgroup', type=int, nargs='+', help="Specify indices of sequences which should not be conserved with the ingroup. Default: none")
    parser.add_argument('--identity', '-i', type=float, default=1, help="Set the minimum identity threshold for conservation. Default: %(default)s")
    parser.add_argument('--properties', '-p', type=float, default=1, help="Set the minimum similarity in chemical properties for conservation analysis. Default: %(default)s")
    parser.add_argument('--gaps', '-g', type=float, default=0.25, help="Set the maximum proportion of gaps allowed in a site for conservation analysis. Default: %(default)s")
    parser.add_argument('--table', '-t', action='store_true', help="Print a summary table instead of a trimmed alignment of the conserved sites.")
    parser.add_argument('--fasta-out', action='store_true', help="Output the trimmed alignment in FASTA format.")
    parser.add_argument('--recode', '-r', choices=recoding_bins.keys(), help="Recode the alignment using the specified scheme before conservation analysis.")
    args = parser.parse_args()
    # Argument validation
    if not os.path.isfile(args.fasta):
        sys.exit("Specified file does not exist: {}".format(args.fasta))
    if (args.identity < 0) or (args.identity > 1):
        sys.exit("Identity threshold must be a proportion between 0 and 1, inclusive.")
    if (args.gaps < 0) or (args.gaps > 1):
        sys.exit("Gaps threshold must be a proportion between 0 and 1, inclusive.")
    if (args.properties < 0) or (args.properties > 1):
        sys.exit("Properties threshold must be a proportion between 0 and 1, inclusive.")
    main(args)

