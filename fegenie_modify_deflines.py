#!/usr/bin/env python3
import sys
import os

def main():
    in_path = sys.argv[1]
    out_path = in_path + ".tmp"
    out = open(out_path, 'w')
    with open(in_path, 'r') as f:
        for line in f:
            out_line = line
            if line[0] == ">":
                sections = line.split("|")
                gene_call = sections[0].lstrip(">")
                contig = sections[1].split(":")[-1].replace("_","")
                out_line = ">" + contig + "_" + gene_call + "\n"
            out.write(out_line)
    out.close()
    os.replace(out_path, in_path)

if __name__ == "__main__":
    main()

