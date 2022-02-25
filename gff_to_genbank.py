#!/usr/bin/env python3
"""Convert a GFF and associated FASTA file into GenBank format.
Modified from https://github.com/chapmanb/bcbb/blob/master/gff/Scripts/gff/gff_to_genbank.py
Usage:
    gff_to_genbank.py <GFF annotation file> <FASTA sequence file>
"""
from __future__ import print_function

import sys
import os

from Bio import SeqIO

# Changed to import from the file GFFParser.py rather than an installed package
import GFFParser

def main(gff_file, fasta_file, molecule_type="DNA"):
    out_file = "%s.gb" % os.path.splitext(gff_file)[0]
    fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    # Modified this line in response to the changed import
    gff_iter = GFFParser.parse(gff_file, fasta_input)
    SeqIO.write(_check_gff(_fix_ncbi_id(gff_iter), molecule_type), out_file, "genbank")
    # Add warning about post-processing
    print("Caution: remember to update relevant fields in the output file.")
    print("For Mauve, add useful info to the SOURCE and ORGANISM fields.")

def _fix_ncbi_id(fasta_iter):
    """GenBank identifiers can only be 16 characters; try to shorten NCBI.
    """
    for rec in fasta_iter:
        if len(rec.name) > 16 and rec.name.find("|") > 0:
            new_id = [x for x in rec.name.split("|") if x][-1]
            print("Warning: shortening NCBI name %s to %s" % (rec.id, new_id))
            rec.id = new_id
            rec.name = new_id
        yield rec

def _check_gff(gff_iterator, molecule_type):
    """Check GFF files before feeding to SeqIO to be sure they have sequences.
    """
    for rec in gff_iterator:
        if "molecule_type" not in rec.annotations:
            rec.annotations["molecule_type"] = molecule_type
        yield _flatten_features(rec)

def _flatten_features(rec):
    """Make sub_features in an input rec flat for output.
    GenBank does not handle nested features, so we want to make
    everything top level.
    """
    out = []
    for f in rec.features:
        cur = [f]
        while len(cur) > 0:
            nextf = []
            for curf in cur:
                out.append(curf)
                if len(curf.sub_features) > 0:
                    nextf.extend(curf.sub_features)
            cur = nextf
    rec.features = out
    return rec

if __name__ == "__main__":
    help_string = """
    Convert a GFF and associated FASTA file into GenBank format.
    Modified from https://github.com/chapmanb/bcbb/blob/master/gff/Scripts/gff/gff_to_genbank.py
    Usage:
        gff_to_genbank.py <GFF annotation file> <FASTA sequence file>
    The output will be written to the same directory as the input GFF file.
    """
    if "-h" in sys.argv or "--help" in sys.argv:
        print(help_string)
        sys.exit(1)
    else:
        main(*sys.argv[1:])
