#!/usr/bin/env python3
import os
import sys
import argparse
import json
import pandas as pd
from collections import defaultdict
from collections import namedtuple

def make_lalign_formatter(df, cols=None):
    """
    Construct formatter dict to left-align columns.

    Parameters
    ----------
    df : pandas.core.frame.DataFrame
        The DataFrame to format
    cols : None or iterable of strings, optional
        The columns of df to left-align. The default, cols=None, will
        left-align all the columns of dtype object

    Returns
    -------
    dict : Formatter dictionary

    """
    if cols is None:
       cols = df.columns[df.dtypes == 'object'] 
    return {col: f'{{:<{df[col].str.len().max()}s}}'.format for col in cols}

def query_cover(hsps, query_length):
    """
    Computes the query coverage of a set of HSPs.
    
    Parameters
    ----------
    hsps : A list containing several HSP dict objects from BLAST JSON output.
    query_length : The integer length of the query sequence.
    
    Returns
    -------
    str : The query coverage formatted as a percent (i.e. '92.4%' rather than 0.9235...)
    """
    # Make a list of all sequence positions in the query that are covered.
    positions = []
    for hsp in hsps:
        positions.extend(range(hsp['query_from'], hsp['query_to'] + 1))
    # Python set objects contain exactly one of each element. Converting our list
    # to a set removes duplicates, and then we can count the covered positions.
    coverage = len(set(positions)) / query_length
    return "{:.1%}".format(coverage)

def print_summary_table(hits, query_len, taxid_flag, length_flag):
    """
    Creates the summary table form of the requested query report.
    
    Parameters
    ----------
    hits : A list containing hit objects from BLAST JSON output.
    query_len : The integer length of the query sequence.
    
    Returns
    -------
    None : Prints output to console
    """
    hit_dict = defaultdict(list)
    # Columns will be: Index, Description, TaxID, Max score, Total score, Query cover, E-value, % ID, Length, Accession
    # Adding to the dictionary in this order means the DataFrame will automatically import the columns in order.
    for hit in hits:
        # There can be multiple HSPs for a given target sequence, meaning there are alignments to
        # different parts of the sequence. The first one will be the best scorer.
        best_hsp = hit['hsps'][0]
        hit_dict['Index'].append(hit['num'])
        # There can be multiple description entries for identical proteins. Just pick the first one.
        hit_dict['Description'].append(hit['description'][0]['title'])
        if taxid_flag:
            hit_dict['TaxID'].append(hit['description'][0]['taxid'])
        # The 'total score' is the bit score for that alignment.
        hit_dict['Bit score'].append(best_hsp['bit_score'])
        # The 'max score' is the sum of all bit scores for all alignments to that protein. 
        hit_dict['Total score'].append(sum([hsp['bit_score'] for hsp in hit['hsps']]))
        hit_dict['Query cover'].append(query_cover(hit['hsps'], query_len))
        hit_dict['E-value'].append(best_hsp['evalue'])
        hit_dict['% Ident.'].append(100*best_hsp['identity']/best_hsp['align_len'])
        if length_flag:
            hit_dict['Align. Len.'].append(best_hsp['align_len'])
        hit_dict['Length'].append(hit['len'])
        hit_dict['Accession'].append(hit['description'][0]['accession'])
    output = pd.DataFrame(hit_dict)
    print(output.to_string(index=False, justify='center', formatters=make_lalign_formatter(output), float_format='{:.4g}'.format))

def split_sequence(hsp, max_width=60):
    """
    Splits the sequence alignment of a HSP into multiple lines,
    replicating the BLAST display format.
    
    Parameters
    ----------
    hsp : An HSP object from the BLAST JSON output
    max_width : An integer specifying the maximum sequence length
                per line. At this point, the alignment will be wrapped.
    
    Returns
    -------
    str : A nicely formatted string representation of the alignment.
    """
    # Split the alignment sequences into fragments of the appropriate length
    split_query = [hsp['qseq'][i:i+max_width] for i in range(0, hsp['align_len'], max_width)]
    split_hit = [hsp['hseq'][i:i+max_width] for i in range(0, hsp['align_len'], max_width)]
    split_midline = [hsp['midline'][i:i+max_width] for i in range(0, hsp['align_len'], max_width)]
    # Format this stuff somehow
    # To get the appropriate placement of the last fragment's end position,
    # it would be easier to make each fragment its own df
    alignment = []
    positions = {'q': hsp['query_from'], 'h': hsp['hit_from']}
    for q, m, h in zip(split_query, split_midline, split_hit):
        data = {}
        data['label'] = ["Query", "", "Sbjct"]
        data['start'] = [positions['q'], "", positions['h']]
        data['seq'] = [q, m, h]
        q_increment = len(q.replace("-","")) - 1
        h_increment = len(h.replace("-","")) - 1
        data['end'] = [positions['q'] + q_increment, "", positions['h'] + h_increment]
        positions['q'] = positions['q'] + q_increment + 1
        positions['h'] = positions['h'] + h_increment + 1
        output = pd.DataFrame(data)
        # The col_space parameter doesn't include inter-column spacing,
        # so this will accommodate sequences up to 9999 positions long.
        alignment.append(output.to_string(index=False, header=False, na_rep="", col_space=4))
    return "\n\n".join(alignment)

def print_alignment_hit(hit):
    """
    Create a nicely formatted report of the alignment(s) between the
    query and the specified hit, loosely replicating BLAST text output.
    
    Parameters
    ----------
    hit : A JSON hit object from the BLAST output.
    
    Returns
    -------
    str : A string representation of the alignment report.
    """
    ## Header for the hit sequence
    # There can be multiple description entries for identical proteins. Just pick the first one.
    print(">{}".format(hit['description'][0]['title']))
    print("Sequence ID: {0}\tLength: {1}".format(hit['description'][0]['id'], hit['len']))
    ## Alignment for each HSP
    for hsp in hit['hsps']:
        print("Range {0[num]}: {0[hit_from]} to {0[hit_to]}\n".format(hsp))
        print("Score: {0[bit_score]:.1f} bits({0[score]}), Expect: {0[evalue]:.3g}".format(hsp))
        stats = {'i': (hsp['identity']/hsp['align_len']), 'p': (hsp['positive']/hsp['align_len']), 'g': (hsp['gaps']/hsp['align_len'])}
        print("Identities: {0[identity]}/{1}({i:.1%}), Positives: {0[positive]}/{1}({p:.1%}), Gaps: {0[gaps]}/{1}({g:.1%})\n".format(hsp, hsp['align_len'], **stats))
        print(split_sequence(hsp))
        print("")

def print_report(report, flags):
    """
    Manage any hit requests and loosely replicate the text output
    from BLAST including headers, tables, and/or alignments.
    
    Parameters
    ----------
    report : A report object from BLAST JSON output.
    flags  : A Namespace object as produced by
             argparse containing input flags.
    
    Returns
    -------
    None : prints output to console
    """
    results = report['results']['search']
    length = results['query_len']
    # Print header information
    print("")
    print("Query title: {0[query_title]}\tQuery ID: {0[query_id]}\tLength: {0[query_len]}".format(results))
    print("Program: {0[version]}\tDatabase: {0[search_target][db]}".format(report))
    if flags.parameters:
        print("Search parameters:")
        for param in report['params']:
            print("{0:>10} : {1}".format(param, report['params'][param]))
    # Select requested hits or all hits if none were specified.
    if flags.hit:
        # We have to do some string casting because json.load will cast to numeric
        # types but command line arguments come in as strings.
        hit_list = [hit for hit in results['hits'] if str(hit['num']) in flags.hit]
        if len(hit_list) != len(flags.hit):
            sys.exit("At least one requested hit is invalid. Please check for typos or use the --help option.")
    else:
        hit_list = results['hits']
    # Print hit information
    if flags.alignments:
        print("Alignments:")
        for hit in hit_list:
            print("")
            print_alignment_hit(hit)
    else:
        print("Sequences producing significant alignments:")
        print_summary_table(hit_list, length, flags.taxid, flags.align_length)
    print("")

def main(args):
    """
    Load the provided JSON results file and handle any report requests.
    If no requests were made, list available reports.
    
    Parameters
    ----------
    args : A Namespace object created by argparse
    
    Returns
    -------
    None : Prints output to console
    """
    ## Load results
    with open(args.input, 'r') as f:
        results = json.load(f)
    # Extract the list of BLAST searches in this results file.
    results = results["BlastOutput2"]
    ## Handle requests from the user
    if args.report is not None: # If the --report option was provided...
        if args.report: # If report titles were specified, print specified reports.
            entries = [entry for entry in results if entry['report']['results']['search']['query_title'] in args.report]
            if len(entries) != len(args.report):
                sys.exit("At least one requested report is invalid. Please check for typos or use the --help option.")
            for entry in entries:
                print_report(entry['report'], args)
        else: # If no report titles were specified, print all reports.
            for entry in results:
                print_report(entry['report'], args)
    else:
        print("Available reports:")
        for entry in results:
            print(entry["report"]["results"]["search"]["query_title"])

if __name__ == "__main__":
    desc = ("A command-line utility to view and navigate downloaded BLAST results files. "
            "Your BLAST results must be downloaded in Single-File JSON format from the NCBI results page."
            "Requires pandas to run.")
    epil = ("Note that if a report title contains spaces, you will need to put it in quotes 'like this' when "
            "passing arguments to the --report option.")
    parser = argparse.ArgumentParser(description=desc, epilog=epil)
    parser.add_argument("-i", "--input", required=True, help="The Single-File JSON output from NCBI BLAST.")
    parser.add_argument("-r", "--report", nargs='*', help="Select one or more reports to view in detail. To view all reports, use this option with no arguments.")
    parser.add_argument("--hit", nargs='+', help="Specify the index of one or more hits to display per report. By default, all hits are displayed.")
    parser.add_argument("-a", "--alignments", action='store_true', help="Display alignment reports instead of a summary table.")
    parser.add_argument("-p", "--parameters", action='store_true', help="Display BLAST input parameters in the header of the results.")
    parser.add_argument("-t", "--taxid", action='store_true', help="Display NCBI taxid in the results summary table.")
    parser.add_argument("-l", "--align-length", action='store_true', help="Display alignment length in the results summary table.")
    args = parser.parse_args()
    if not os.path.isfile(args.input):
        sys.exit("Specified file does not exist: {}".format(args.input))
    main(args)

