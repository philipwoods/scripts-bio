#!/usr/bin/env python3
## Import required modules
import sys
import os
import xml.etree.ElementTree as ET

def main():
    ## Global variables
    help_string = """
    Synopsis:
        python refseq_exclusion.py [-h] [--version] [<output options>] input-file

    Description:
        This utility extracts which assemblies in downloaded NCBI data were excluded 
        from RefSeq and the reason they were excluded. You must provide the path to
        an assembly_result.xml file corresponding to your NCBI data. The XML file can
        be downloaded from the 'Send to:' menu in the upper right of your NCBI search
        results page with the 'File' destination and 'XML' format options selected.
    
    Options:
        -h | --help
            Display this help text and exit.
        --version
            Display version information and exit.
        -l | --list <identifier>
            Outputs a simple list of the specified identifier of excluded assemblies.
        -r | --reason <phrase>
            Filters output to only show assemblies excluded for the specified reason.
    """
    version_string = """
    Last updated June 21 2019
    """
    flags = {
        'list': False,
        'reason': False
    }
    options = {
        'identifier': "",
        'reason': ""
    }

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
        elif arg in ["-l", "--list"]:
            flags['list'] = True
            options['identifier'] = args.pop(0)
        elif arg in ["-r", "--reason"]:
            flags['reason'] = True
            options['reason'] = args.pop(0)
        elif arg not in ["-", "--"]:
            sys.exit("Invalid option: " + arg)

    if len(args) < 1:
        sys.exit("Missing input file argument.\nUse the --help option to learn more.")
    if len(args) > 1:
        sys.exit("Extra argument(s) detected.\nUse the --help option to learn more.")

    target = args.pop(0)
    if not os.path.isfile(target):
        sys.exit("Specified file does not exist: " + target)

    if not flags['list']:
        print("Genome assemblies excluded from RefSeq")
        print("Target file: " + os.path.abspath(target))
        if flags['reason']: print("Specified reason: " + options['reason'])
        print("------------------------------------------------")

    find_exclusions(target, flags, options)

def find_exclusions(target_file, flags, arguments):
    with open(target_file,'r') as f:
        # Add a <data> tag as a wrapper to provide a single root node.
        text = "<data>" + f.read() + "</data>"
    root = ET.fromstring(text)
    del text # So we don't keep two copies of the file in memory
    
    ## Iterate through files looking for excluded assemblies.
    counter = 0
    for DS in root.findall("DocumentSummary"):  # For each DocumentSummary node...
        exclusion = DS.find("ExclFromRefSeq")  # Find the exclusion node
        if exclusion.text is not None:  # If the node has any contents...
            # Get a list of the exclusion reasons
            reason_list = [reason.text for reason in exclusion.findall("string")]
            # If we're filtering and our specified reason isn't in the list, skip
            if flags['reason'] and (arguments['reason'] not in reason_list): continue
            # Format output depending on the flags
            if flags['list']:
                print(DS.find(arguments['identifier']).text)
            else:
                counter += 1
                print("Organism: " + DS.find("Organism").text)
                print("Assembly Name: " + DS.find("AssemblyName").text)
                print("Accession: " + DS.find("AssemblyAccession").text)
                print("Reasons for exclusion from RefSeq:")
                for reason in reason_list:
                    print("    " + reason)
                print("")
    if not flags['list']: print("Total excluded entries: " + str(counter))

if __name__ == "__main__":
    main()
