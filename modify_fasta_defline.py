#!/usr/bin/env python3
## Import required modules
import os
import sys

def main():
    help_string = """
    Synopsis:
        python modify_fasta_defline.py FASTA_FILE GENOME_NAME INFO_FILE
    
    Description:
        This is an utility for the script anvi-script-get-sequences-by-function.sh
        and is not generally intended for independent use. For that script to work,
        both scripts must be in the same directory. Assumes exactly one record in
        the info file per gene ID.

        Lines in the INFO_FILE should be formatted as ID|info1|info2 etc.
    """
    version_string = """
    Last updated 14 September 2020
    """

    ## Handle arguments
    if '-h' in sys.argv or '--help' in sys.argv:
        print(help_string)
        sys.exit(1)
    if '-v' in sys.argv or '--version' in sys.argv:
        print(version_string)
        sys.exit(1)
    fasta_path = sys.argv[1]
    genome_name = sys.argv[2]
    info_path = sys.argv[3]
    temp_path = fasta_path + ".tmp"

    out_list = []

    with open(info_path) as f:
        info_list = f.readlines()
    info_dict = {}
    for line in info_list:
        gene_id = line.partition("|")[0]
        info_dict[gene_id] = line

    print("Adding information to FASTA defline(s)...")
    with open(fasta_path) as f:
        for line in f:
            if line[0] == ">":
                gene_id = line[1:].rstrip()
                # Output NO ANNOTATION if the gene ID isn't in the function file
                gene_info = info_dict.get(gene_id)
                if gene_info:
                    out_line = ">" + genome_name + "|" + gene_info
                else:
                    out_line = ">" + genome_name + "|" + gene_id + "|NO ANNOTATION\n"
                out_list.append(out_line)
            else:
                out_list.append(line)

    print("Cleaning up...")
    out_string = "".join(out_list)
    with open(temp_path, 'w') as f:
        f.write(out_string)
    os.replace(temp_path, fasta_path)

if __name__ == "__main__":
    main()
