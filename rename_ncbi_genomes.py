#!/usr/bin/env python3
## Import required modules
import sys
import os
import xml.etree.ElementTree as ET

def main():
    ## Global variables
    help_string = """
    Synopsis:
        python rename_genomes.py [-h] [--version] [-v] database xml-file assembly-dir

    Description:
        This utility helps give NCBI genome assembly files more meaningful names. You
        must specify GenBank or RefSeq as the source database. GenBank is the default.
        For brevity, you may also use the first letter of each database name to do
        this. You must also provide the path to an XML file formatted like a NCBI
        assemly_result.xml file. Finally, you must provide the path to a directory
        containing .fna files named using the default NCBI naming scheme.

        The XML file can be downloaded from the 'Send to:' menu in the upper right of
        your NCBI search results page with the 'File' destination and 'XML' format
        options selected.
    
    Options:
        -h | --help
            Display this help text and exit.
        --version
            Display version information and exit.
        -v | --verbose
            Provide more extensive printed output.
    """
    version_string = """
    Last updated June 21 2019
    """
    verbose_flag = False
    database = ""

    ## Handle arguments
    args = sys.argv[1:]
    while len(args) > 0 and args[0][0] == "-":
        arg = args.pop(0)
        if arg == "-h" or arg == "--help":
            print(help_string)
            sys.exit(0)
        elif arg == "--version":
            print(version_string)
            sys.exit(0)
        elif arg in ["-v", "--verbose"]:
            verbose_flag = True
        elif arg not in ["-", "--"]:
            sys.exit("Invalid option: " + arg)

    if len(args) < 3:
        sys.exit("Missing argument.\nUse the --help option to learn more.")
    if len(args) > 3:
        sys.exit("Extra argument detected.\nUse the --help option to learn more.")

    db_arg = args.pop(0)
    if db_arg in ["g", "G", "GenBank"]:
        database = "Genbank"
    elif db_arg in ["r", "R", "RefSeq"]:
        database = "RefSeq"
    else:
        err_str1 = "You must specify GenBank or RefSeq as the source genome database.\n"
        err_str2 = "Use the --help option to learn more."
        err_str = err_str1 + err_str2
        sys.exit(err_str)

    xml_file = args.pop(0)
    if not os.path.isfile(xml_file):
        sys.exit("Specified file does not exist: " + xml_file)

    target_dir = args.pop(0)
    if not os.path.isdir(target_dir):
        sys.exit("Specified directory does not exist: " + target_dir)

    if verbose_flag:
        print("Current directory: " + os.getcwd())
        print("Selected database: " + database)
        print("Target XML file: " + xml_file)
        print("Target directory: " + target_dir)
    
    rename_files(xml_file, target_dir, database, verbose_flag)

def rename_files(xml_file, directory, db, verbose):
    ## Input validation and parsing
    directory_file_list = list(filter(lambda x: x[-4:]==".fna", os.listdir(directory)))
    if verbose:
        print("".join(["Found ", str(len(directory_file_list)), " assembly files in the target directory."]))
        print("Reading XML data...")
    with open(xml_file,'r') as f:
        # Add a <data> tag as a wrapper to provide a single root node.
        text = "<data>" + f.read() + "</data>"
    root = ET.fromstring(text)
    del text # So we don't have two copies of the file in memory
    if verbose: print("Found " + str(len(root)) + " entries in assembly_result.xml")
    if len(root) != len(directory_file_list):
        sys.exit("Mismatched number of files and file descriptions.")
    
    ## Iterate through and rename files
    for DS in root.findall("DocumentSummary"): # For each DocumentSummary node in the file...
        synonym = DS.find("Synonym").find(db).text
        assembly_name = DS.find("AssemblyName").text
        species = DS.find("SpeciesName").text
        accession = DS.find("AssemblyAccession").text
    
        old_name = "_".join([synonym, assembly_name, "genomic.fna"])
        old_name = "".join(x if x not in """:\\/<'">%?|* """ else "_" for x in old_name)
        if not os.path.isfile(directory + old_name):
            err_str1 = "Requested file does not exist: " + old_name + "\n"
            err_str2 = "You may have already renamed these files."
            err_str = err_str1 + err_str2
            sys.exit(err_str)
        new_name = "".join([species, "_", accession, ".fna"])
        new_name = "".join(x if x not in """:\\/<'">%?|* """ else "_" for x in new_name)
        if verbose: print("".join(["Renaming ", old_name, " to ", new_name]))
        os.rename(directory + old_name, directory + new_name)
    print("Done")

if __name__ == "__main__":
    main()
