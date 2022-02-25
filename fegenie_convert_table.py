#!/usr/bin/env python3
import sys
import os.path
import csv

def main():
    in_path = sys.argv[1]
    directory = os.path.dirname(in_path)
    out_path = os.path.join(directory, "fegenie-anvi-table.tsv")
    out = open(out_path, 'w')
    print("gene_callers_id\tsource\taccession\tfunction\te_value", file=out)
    with open(in_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_caller_id = row['orf'].split("_")[-1]
            source = "FeGenie"
            accession = "NoAccession"
            category = row['category'].capitalize().replace("_"," ").replace("-",": ")
            hmm = row['HMM'].replace("-"," ")
            function = ", ".join((category, hmm))
            e_value = "0"
            print("\t".join((gene_caller_id, source, accession, function, e_value)), file=out)
    out.close()

if __name__ == "__main__":
    main()

