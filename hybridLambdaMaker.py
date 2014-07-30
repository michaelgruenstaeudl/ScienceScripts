#!/usr/bin/env python2
'''Converting a left-ladderized phylogenetic tree into input for Hybrid-Lambda'''
__author__ = "Michael Gruenstaeudl, PhD"
__copyright__ = "Copyright (C) 2014 Michael Gruenstaeudl"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.07.29.1300"
__status__ = "Pseudocode"

#########################
### IMPORT OPERATIONS ###
#########################

from termcolor import colored
import argparse, os, sys
import GeneralStringOperations as GSO

###################
### DEFINITIONS ###
###################

def addNodeNumbers(intree):
    ''' Adding node numbers to tree.
    Args:
        intree: a left-ladderized tree string without node labels.
    Returns:
        outtree: a left-ladderized tree string with node labels.
    '''
    outtree = []
    nNodes = (intree.count(")")*2)+2
    for c in intree:
        if c == ")":
            outtree.append(c + "s" + str(nNodes))
            nNodes -= 1
        else:
            outtree.append(c)
    return ''.join(outtree)


def main(treeName, parentInfo):
    # Reading plain tree
    tree = open(treeName, "r").read()

    # Parsing parentInfo into dictionary
    aDict = {}
    for i in parentInfo.split(";"):
        aDict[i.split(",")[0]] = i.split(",")[1]
    # Hypothetical content of aDict at this point:
    #    {"B":"0.6", "C":"0.4"}

    # Adding hybrid information to intree
    addHybrInfo(intree, hybInfoDict)
#    tree = addHybrInfo(intree, hybInfoDict)

    # Adding nodes to tree
#    treeNN = addNodeNumbers(tree)

    # Saving output to file
#    outf = open(GSO.rmext(pathAndInfilename)+".hybLambda.INPUT", "w")
#    outf.write(''.join(outd))


def addHybrInfo(intree, hybInfoDict):
    ''' Adds hybrid information to tree string in accordance with the 
        requirements of hybrid-Lambda.
    Args:
        intree: a tree string without hybrid info.
        hybInfoDict: dictionary with hybrid parent info of format {"B":"0.6", "C":"0.4"}

    Returns:
        outtree: a tree string with hybrid info.
    '''

    outtree = []
    # Looping over the input taxa and incorporating the hyrbid info for each 
    # taxon in list
    for key,val in hybInfoDict.iteritems():

        # 1. Replace parent's brlen with adjusted brlen
        if not intree.find(key):
            print colored("  Error: Specified hybrid parent not found in \
input tree.", 'red')
            break()
        else:
            # Parse out the branch length immediately following the key, 
            # save as "brlen"
            brlen = float(GSO.exstr(intree,key+":",")"))
            # Replace said branch length with half of value
            intree = GSO.replstr(intree,key+":",")",str(brlen*0.5))

        # 2. Add hybrid info to parent
        # Split intree into three sections
        split1 = GSO.csplit(intree,key,rightflag=T)
        aList = [split1[0], GSO.csplit(split1[1],key+":"+str(brlen*0.5),rightflag=F)
        print aList

    ## add [open parenthesis "(" + "h#" + hybr. prob. + colon + 1/2brlen + comma] to infstr.find(c)-1

    ## add closed parenthesis ")" to infstr.afind(")",c)


###############
### EXECUTE ###
###############

print ""
print colored("  Script name: "+sys.argv[0], 'cyan')
print colored("  Author: "+__author__, 'cyan')
print colored("  Version: "+__version__, 'cyan')
print ""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converting a left-ladderized phylogenetic tree into input for Hybrid-Lambda; '+__copyright__)
    parser.add_argument('-t','--tree', help='name of input tree (left-ladderized!)', default="infile.tre", required=True)
    parser.add_argument('-p','--parents', help='info on parental taxa; format: <parent1>,<likelih.parent1>;<parent2>,<likelih.parent2>', default="A,0.6;B,0.4", required=True)
    args = parser.parse_args()

main(args.tree, args.parents)

print ""
print colored("  Done.", 'cyan')
print ""

