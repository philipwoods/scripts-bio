#!/usr/bin/env python3
## Import required modules
import sys
import os
import string

def main():
    ## Global variables
    help_string = """
    Synopsis:
        python gen-anvi-external-genomes.py [-h] [--version] [-o <out-path>] [-d <dir>]
                                                                             [-f <file>]
                                                                             [<db1> ...]

    Description:
        This script creates an anvi'o external genomes file for the specified contig
        databases. This external genomes file can be used as an argument for anvi'o
        tools with the --external-genomes option. The databases can be specified as
        a list of file names, a directory containing .db files, or a file listing
        the names of the database files of interest.

    Options:
        -h | --help
            Display this help text and exit.
        --version
            Display version information and exit.
        -o | --output-path <output-path>
            Specify a path where an output file should be written.
        -d | --directory <directory>
            Specify a directory for the script to look in for database files.
        -f | --file <file>
            Specify a file containing the names of the desired database files.
    """
    version_string = """
    Last updated July 8 2019
    """
    flags = {
        'file_out': False,
        'input': "list"
    }
    options = {
        'in_files': None,
        'in_path': None,
        'out_path': None
    }

    ## Handle arguments
    args = sys.argv[1:]
    while len(args) > 0 and args[0][0] == "-":
        arg = args.pop(0)
        if arg in ["-", "--"]:
            break
        elif arg == "-h" or arg == "--help":
            print(help_string)
            sys.exit(0)
        elif arg == "--version":
            print(version_string)
            sys.exit(0)
        elif arg in ["-d", "--directory"]:
            flags['input'] = "dir"
            options['in_path'] = args.pop(0)
        elif arg in ["-f", "--file"]:
            flags['input'] = "file"
            options['in_path'] = args.pop(0)
        elif arg in ["-o", "--output-path"]:
            flags['file_out'] = True
            options['out_path'] = args.pop(0)

    # If the input is a list of files (default), grab that list.
    if flags['input'] == "list":
        options['in_files'] = args
    # If the input is anything else but there are still arguments, throw an error
    elif len(args) > 0:
        sys.exit("Extra argument(s) detected.\nUse the --help option to learn more.")
    # If the input is a directory, get the file list
    if flags['input'] == "dir":
        options['in_files'] = filter(lambda x: os.path.splitext(x)[1]==".db", os.listdir(options['in_path']))
        # os.listdir() gives a list of files as if you were in the directory, not paths to them
        options['in_files'] = [os.path.join(options['in_path'], x) for x in options['in_files']]
    # If the input is a file, get the file list
    if flags['input'] == "file":
        with open(options['in_path'], 'r') as f:
            options['in_files'] = [x.rstrip("\n") for x in f.readlines()]

    # Checking for errors
    if flags['input'] == "file" and not os.path.isfile(options['in_path']):
        sys.exit("Specified file does not exist: " + options['in_path'])
    if flags['input'] == "dir" and not os.path.isdir(options['in_path']):
        sys.exit("Specified directory does not exist: " + options['in_path'])

    build_file(flags, options)

def build_file(flags, arguments):
    # Select the appropriate output based on user arguments
    out_file = sys.stdout
    if flags['file_out']:
        out_file = open(arguments['out_path'], 'w')

    # Write header to file
    out_file.write("name\tcontigs_db_path\n")
    # The file needs two columns: name and contigs_db_path. The name is a handle for
    # anvi'o to use for a genome, and the contigs_db_path is the path to the contigs
    # database file for that genome. The Anvi'o documentation recommends only using
    # alphanumerics and underscores in names.
    allowed_name_chars = "".join([string.ascii_letters, string.digits, "_"])
    # Iterate over file list and build output file
    for filepath in arguments['in_files']:
        if not os.path.splitext(filepath)[1] == ".db":
            sys.exit("Filename must end in .db: " + filepath)
        # Get the file name (after the last pathname separator) excluding the extension
        filename_root = os.path.splitext(os.path.split(filepath)[1])[0]
        # Replace everything disallowed with an underscore
        name = "".join(x if x in allowed_name_chars else "_" for x in filename_root)
        # Output the absolute path of the file
        contigs_db_path = os.path.abspath(filepath)
        out_file.write(name + "\t" + contigs_db_path + "\n")
    out_file.close()

if __name__ == "__main__":
    main()
