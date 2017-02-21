#!/usr/bin/env python2
'''PROTOCOL 001 - script 01'''
__author__ = "Michael Gruenstaeudl, PhD"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.09.07.2100"
__status__ = "Working"


#####################
# IMPORT OPERATIONS #
#####################


import argparse
import dendropy
import GeneralStringOperations as GSO
import sys
from termcolor import colored


########
# MAIN #
########


def main(msCmd, inTree):

# Get pop info from commandline
    # read commands
    handle = open(msCmd).read()
    # extract relevant section
    popLst = GSO.exstr(handle, "-I", "-")
    # remove leading and trailing whitespaces, make into list
    popLst = popLst.lstrip().rstrip().split()
    # remove:   first element (merely number of pops, not size)
    #           last element: the outgroup
    popLst = popLst[1:-1]

# Generate dictionary for label replacement
    lttrs = map(chr, range(97, 123))
    # remove "o" from list
    lttrs.remove("o")
    # add a zero element to the beginning of lttrs
    lttrs = ["0"] + lttrs
    # Initializations
    Sum = 0
    replaceDict = {}
    # loop over popLst
    for cntr, pop in enumerate(popLst, start=1):
        for num in range(Sum+1, Sum+int(pop)+1):
            replaceDict[str(num)] = lttrs[cntr] + str(num).zfill(4)
        Sum += int(pop)
    # include an outgroup
        replaceDict[str(Sum+1)] = "o" + str(Sum+1).zfill(4)

# Load trees and perform label replacement
    trees = dendropy.TreeList.get_from_path(inTree, 'newick')
    # leaf_nodes() need to be changed only for the first tree of trees;
    # the taxonnames of the other tree of trees are changed simultaneoulsy
    for leaf in trees[0].leaf_nodes():
        leaf.taxon.label = replaceDict[leaf.taxon.label]

# Save trees back to file
    trees.write_to_path(GSO.rmext(inTree)+".relbld.tre", 'newick')


###########
# EXECUTE #
###########

print ""
print colored("    Script name: "+sys.argv[0], 'cyan')
print colored("    Author: "+__author__, 'cyan')
print colored("    Version: "+__version__, 'cyan')
print ""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='PROTOCOL 001 - script 01')

    parser.add_argument('-m',
                        '--mscmd',
                        help='name of mscmd file',
                        default="msCmd.txt",
                        required=True)

    parser.add_argument('-t',
                        '--tree',
                        help='name of treefile',
                        default="intree.tre",
                        required=True)

    args = parser.parse_args()

main(args.mscmd, args.tree)

print colored("    Done.", 'cyan')
print ""
