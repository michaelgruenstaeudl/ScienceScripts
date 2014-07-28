#!/usr/bin/env python2
'''Parsing a phylogenetic tree file to be compatible as input to Hybrid-Lambda'''
__author__ = "Michael Gruenstaeudl, PhD"
__copyright__ = "Copyright (C) 2014 Michael Gruenstaeudl"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.07.28.1400"
__status__ = "Pseudocode"

#########################
### IMPORT OPERATIONS ###
#########################

import os
import sys

print("I require left-ladderizied trees.")

if len(sys.argv) < 2:
    print("python hybridLambdaPrep.py inputfile")
    sys.exit()

inf = open(sys.argv[1], "r")
infstr = inf.read()
outf = open(sys.argv[1]+".out", "w")

# Adding node numbers to tree
outd = []
count = 2*infstr.count(")")+1
for c in infstr:
    if c == ")":
        outd.append(c + "s" + str(count))
        count -= 1
    else:
        outd.append(c)

# Looping over the input taxa and incorporating the hyrbid info for each taxon in list
for taxon in adict.iteritems():
    out = AddHybInf(infstr, taxon[0], taxon[1])

# Hypothetical content of adict:
#    {"A":"0.9", "B":"0.1"}


def addHybrInfo(infstr, taxonName, taxonName):
    """ Adds hybrid information to tree string in accordance with the 
        requirements of hybrid-Lambda.

    Args:
        infstr: A string without hybrid info.
        taxonName: name of putative hybrid parent.
        taxonName: probabiity of hybridisation.

    Returns:
        outstring: A string containing the new hybrid info.
    """

# Adding hybrid specs to tree
# For c in infstr:
#   if c = keyw:
#       parse out the branch length immediately following the keyw
#       add [open parenthesis + "h#" + hybr. prob. + colon + brlen + comma] to infstr.find(c)-1
#       add [closed parenthesis] to infstr.afind(")",c)


outf.write(''.join(outd))

