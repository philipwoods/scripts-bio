#!/usr/bin/env python3
## Import required modules
import sys
import os
import xml.etree.ElementTree as ET

def main():
    ## Global variables
    help_string = """
    Synopsis:
        python get_ncbi_names.py [-h] [--version] xml-file

    Description:
        This utility extracts the filename that would correspond to a given NCBI
        genome file from an XML file formatted like a NCBI assemly_result.xml file.

        The XML file can be downloaded from the 'Send to:' menu in the upper right of
        your NCBI search results page with the 'File' destination and 'XML' format
        options selected.
    
    Options:
        -h | --help
            Display this help text and exit.
        --version
            Display version information and exit.
    """
    version_string = """
    Last updated June 21 2019
    """

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
        elif arg not in ["-", "--"]:
            sys.exit("Invalid option: " + arg)

    if len(args) > 1:
        sys.exit("Extra argument detected.\nUse the --help option to learn more.")

    xml_file = args.pop(0)
    if not os.path.isfile(xml_file):
        sys.exit("Specified file does not exist: " + xml_file)

    list_files(xml_file)

def list_files(xml_file):
    ## Input validation and parsing
    with open(xml_file,'r') as f:
        # Add a <data> tag as a wrapper to provide a single root node.
        text = "<data>" + f.read() + "</data>"
    root = ET.fromstring(text)
    del text # So we don't have two copies of the file in memory
    
    ## Iterate through and list file names
    for DS in root.findall("DocumentSummary"): # For each DocumentSummary node in the file...
        species = DS.find("SpeciesName").text
        accession = DS.find("AssemblyAccession").text
        new_name = "_".join([species, accession])
        new_name = "".join(x if x not in """:\\/<'">%?|* """ else "_" for x in new_name)
        print(new_name)

if __name__ == "__main__":
    main()
