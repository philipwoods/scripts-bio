#!/usr/bin/env python
import os
import sys
import argparse

def main():
    # Help and setup information
    sequence_file_types = [".fa", ".fasta", ".fna"]
    annotation_file_types = [".gff", ".gff3"]
    func_desc = """
    Create definitions for one or more Artemis projects as formatted in the .artemis.project.properties
    configuration file. If one or more directories are provided, creates one project description per
    directory using their contents. If a file is provided, each line of the file will be used to create
    a separate project description.

    Basic Artemis projects ask for a project name / title, a nucleotide sequence, and one or more
    annotations. The title and project name will be taken from the directory name. The rest will be
    assumed based on file type.
    """
    type_list = "Sequence file extensions:\n" + "\n".join(sequence_file_types) + "\n\nAnnotation file extensions:\n" + "\n".join(annotation_file_types)

    # Handle arguments
    parser = argparse.ArgumentParser(description=func_desc, epilog=type_list, allow_abbrev=False, formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-d", "--dir", nargs='+',
        help="Directory containing files for an Artemis project. Multiple can be specified.")
    group.add_argument("-f", "--file", type=argparse.FileType(),
        help="File listing directories for Artemis projects.")
    parser.add_argument("-o", "--out", nargs='?', type=argparse.FileType('w'), default=sys.stdout,
        help="Specify output file. If not provided, will write to stdout.")
    args = parser.parse_args()

    # Consolidate different input sources
    projects = []
    if args.file:
        projects = args.file.readlines()
    else:
        projects = args.dir

    # Create project entries
    for project in projects:   
        if not os.path.isdir(project):
            sys.exit("Directory does not exist: {}".format(project))
        # Get consistent path format, then take the directory name
        name = os.path.basename(os.path.abspath(project))
        # Get the absolute paths of the sequence and annotation files
        files = [os.path.abspath(f) for f in os.listdir(project) if os.path.isfile(f)]
        seqs = [f for f in files if (os.path.splitext(f) in sequence_file_types)]
        annotations = [f for f in files if (os.path.splitext(f) in annotation_file_types)]
        args.out.write("#\n")
        args.out.write("project.{0}.title={0}\n".format(name))
        for seq in seqs:
            args.out.write("project.{0}.sequence={1}".format(name, seq))
        for annotation in annotations:
            args.out.write("project.{0}.annotation={1}".format(name, annotation))

if __name__ == "__main__":
    main()

