#!/usr/bin/env python3
import os
import sys
import argparse
import textwrap
import pandas as pd

def main(args):
    # Handle input
    out_base = os.path.abspath(args.output)
    genes_df = pd.read_csv(args.genes, sep=args.separator, dtype=str)
    genes = genes_df.columns[2:] # Exclude 'name' and 'contigs_db_path'
    out_files = [os.path.join(out_base, gene + ".fa") for gene in genes]
    for f in out_files:
        if args.force: # Erase any existing files with these names
            blank = open(f, 'w')
            blank.close()
        elif os.path.exists(f):
            sys.exit("File already exists: {}\nIf you want to overwrite, use the --force option.".format(f))
    # Begin exporting the requested genes
    print("Exporting genes...")
    for row in genes_df.itertuples():
        for index, gene in enumerate(row._fields):
            # Each row namedtuple will have the fields 'Index', 'name', 'contigs_db_path', followed by the gene names
            # We don't want to iterate over anything but the gene names
            if gene in ['Index', 'name', 'contigs_db_path']:
                continue
            if pd.isna(row[index]): # If we aren't asking for this gene from this contigs db...
                continue
            # Set up output paths
            sequences = os.path.join(out_base, gene + ".fa")
            temp = os.path.join(out_base, gene + "_" + row.name + ".fa")
            os.system("anvi-get-sequences-for-gene-calls -c {contigs} --gene-caller-ids {id} --get-aa-sequences -o {out}".format(contigs=row.contigs_db_path, id=row[index], out=temp))
            # Append the individual sequence to the gene file.
            print("Adding sequence to {}".format(sequences))
            with open(sequences, 'a') as f, open(temp, 'r') as t:
                contents = t.read()
                # Existing defline should be >GENE_CALLERS_ID
                defline = ">{organism} gene:{gene_name}|gene_callers_id:".format(organism=row.name, gene_name=gene)
                contents = contents.replace(">", defline)
                f.write(contents)
            print("Removing temp file {}".format(temp))
            os.remove(temp)
    print("Done exporting genes.")
    # Begin aligning the collected gene files
    print("Aligning genes...")
    for f in out_files:
        aligned = os.path.splitext(f)[0] + "-aligned.fa"
        os.system("muscle -in {input} -out {output} -quiet".format(input=f,output=aligned))
    print("Done aligning genes.")

if __name__ == "__main__":
    # Handle arguments
    desc = """
    Takes in a tab-delimited file specifying genomes and genes export from each and align.
    The gene names must be unique, and all names should not contain whitespace. If there are
    multiple databases corresponding to the same genome, you can provide the same name in
    multiple rows. Some gene ID fields can be left empty.
    
    You must have an anvio environment active to run this.
    """
    epi = """
    The input file must be formatted as follows:
    
      name        contigs_db_path <gene_name_1>   <gene_name_2>   . .
      <name_1>    <path_1>        <gene_id_1,1>   <gene_id_1,2>   . .
      <name_2>    <path_2>        <gene_id_2,1>   <gene_id_2,2>
      .                                                           .
      .                                                             .
    
    Note that the first two columns are identical to an anvio external genomes file.
    """
    parser = argparse.ArgumentParser(description=textwrap.dedent(desc), epilog=textwrap.dedent(epi), formatter_class=argparse.RawDescriptionHelpFormatter)
    #parser.add_argument('-c', '--concatenate', action='store_true', help="Create a concatenated alignment of the genes. Default is off.")
    parser.add_argument('-f', '--force', action='store_true', help="Overwrite existing files. Default is off.")
    parser.add_argument('-s', '--separator', default='\t', help="The separator used in the input genes file. Default is TAB.")
    parser.add_argument('-g', '--genes', required=True, help="Gene definition file. See epilog for format information.") 
    parser.add_argument('-o', '--output', required=True, help="Specify output directory.")
    args = parser.parse_args()
    # Validate input
    if not os.path.isfile(args.genes):
        sys.exit("The specified file does not exist: {}".format(os.path.abspath(args.genes)))
    # Check that we have the appropriate environment active
    if 'anvio' not in os.getenv('CONDA_DEFAULT_ENV'):
        sys.exit("You must activate an anvio environment before running this script.")
    # Run main function
    main(args)

