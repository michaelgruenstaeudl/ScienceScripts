#!/usr/bin/env python2
''' GenBank data formatting '''
__author__ = "Michael Gruenstaeudl, PhD"
__copyright__ = "Copyright (C) 2014 Michael Gruenstaeudl"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.12.13.0100"
__status__ = "Working"
__example__ = "python2 this_script.py -f sequences.fas -t infotable.tsv"

# Abbreviations:
# Source Modifier = sm

# IMPORT OPERATIONS
from termcolor import colored
import argparse
import sys


def main(fasta_filename, tsv_filename):
    indict, outdict = {}, {}

    tsv_file = open(tsv_filename, "r").readlines()
    for line in tsv_file[1:]:
        tmp_list = line.split("\t")
        indict[tmp_list[0]] = tmp_list
        for key in indict:
            indict[key] = [e.rstrip() for e in indict[key]]                
            sm_seqid = indict[key][0]
            sm_organism = " [organism=" + indict[key][1] + "]"
            sm_authority = " [authority=" + indict[key][2] + "]"
            sm_country = " [country=" + indict[key][4] + "]"
            sm_specvoucher = " [specimen-voucher=" + indict[key][6] + ": " + indict[key][5] + "]"
            sm_variety = " [variety=" + indict[key][7] + "]"

            outdict[key] = [sm_seqid, sm_organism, sm_authority,
                            sm_country, sm_specvoucher, sm_variety]

    fasta_file = open(fasta_filename, "r").read()
    for key in outdict:
        outstring = ''.join(outdict[key])
        outstring = outstring.replace("[variety=]","").replace("[note=]","").rstrip()
        fasta_file = fasta_file.replace(key, outstring)

        if fasta_file.find(key) == -1:
            print colored("  Not included in outfile: ", 'magenta'), key
            
    outfile = open("OUT_" + fasta_filename[:-4] + ".fsa", "w")
    outfile.write(fasta_file)           
    outfile.close()


# EXECUTE

print ""
print colored("  Script name: " + sys.argv[0], 'cyan')
print colored("  Author: " + __author__, 'cyan')
print colored("  Version: " + __version__, 'cyan')
print ""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GenBank data formatting; \
                                     2014 Michael Gruenstaeudl')
    parser.add_argument('-f', '--fasta_file',
                        help='/path_to_working_dir/sequences.fas',
                        required=True)
    parser.add_argument('-t', '--tsv_file',
                        help='/path_to_working_dir/infotable.tsv',
                        required=True)
    args = parser.parse_args()

main(args.fasta_file, args.tsv_file)

print ""
print colored("  Done.", 'cyan')
print ""


