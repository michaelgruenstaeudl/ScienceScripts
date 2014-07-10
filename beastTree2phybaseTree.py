#!/usr/bin/env python2
'''Converting a species tree string from the *BEAST format to the PHYBASE format.'''
__author__ = "Michael Gruenstaeudl, PhD"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.07.09.2300"
__status__ = "Working"

#########################
### IMPORT OPERATIONS ###
#########################

from termcolor import colored
import argparse, dendropy, numpy, sys
sys.path.insert(0,"/home/michael/git/ScienceScripts/")
import GeneralStringOperations as GSO

########################
### GLOBAL VARIABLES ###
########################
# none

###################
### DEFINITIONS ###
###################

def changing_location_of_metadata(pathAndInfilename):
    ''' Loading a set of trees in *BEAST format and saving them to
        phybase format. Simultaneously converting dmv values to theta value. '''
    # Loading *BEAST file into DendroPy
    inData = dendropy.DataSet.get_from_path(pathAndInfilename, \
    "beast-summary-tree", extract_comment_metadata=True)
    # Looping over all trees
    for tree in inData.tree_lists[0]:
        nodes = tree.nodes()
        # Looping over every node of a given tree
        for node in nodes:
            handle = node.annotations.get_value("dmv")
            # Decision tree regarding different BEAUTi versions
            # For BEAUTi v1.7 and lower, b/c dmv = '0.15':
            if type(handle) is str:
                theta = handle
            # For BEAUTi v1.8 and higher, b/c dmv = ['0.1','0.2']:
            if type(handle) is list:
                theta = converting_dmv_to_theta(handle)
            # Dropping existing metadata
            node.annotations.drop()
            # Saving new metadata as replacement for dropped metadata
            node.annotations.add_new("theta",theta)
    # Returning tree as string
    outData = inData.as_string("nexus")
    return outData

def converting_dmv_to_theta(inList):
    ''' Converting the dmv values (always 2 per node) to a theta value
        inList: dmv = ['0.1','0.2'] '''
    # Converting strings to floats in list
    tmpList = [float(number) for number in inList]
    # Calculating artithmetic mean
    tmpValue = numpy.mean(tmpList)
    # Multiplying by two
    outValue = str(float(tmpValue)*2)
    return outValue

def formatting_tree_string(inString):
    ''' Formatting a tree string '''
    # Convert entire string into lower case
    handle = inString.lower()
    # Remove "taxa" section of nexus file
    handle = handle.replace(handle[handle.find("begin taxa;"):handle.find("end;", \
    handle.find("begin taxa;"))+4],"")
    handle = handle.replace("[&theta=","#")
    handle = handle.replace("[&r]","")
    # Note: Remove "[&r]" prior to "]"; otherwise "[&r]" is no longer found.
    handle = handle.replace("]","")
    return handle

############
### MAIN ###
############

def main(pathAndInfilename):
    # Conducting main functions
    treeList = changing_location_of_metadata(pathAndInfilename)
    outData = formatting_tree_string(treeList)
    # Writing outfile
    outfileName = GSO.rmext(pathAndInfilename)+".phyb"
    outFile = open(outfileName,"w")
    outFile.write(outData)
    outFile.close()

###############
### EXECUTE ###
###############

print colored("  Script name: "+sys.argv[0], 'cyan')
print colored("  Author: "+__author__, 'cyan')
print colored("  Version: "+__version__, 'cyan')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converting a species tree string from the *BEAST format to the PHYBASE format; 2014 Michael Gruenstaeudl')
    parser.add_argument('-i','--pathAndInfilename', help='/path_to_working_dir/infilename', default="~/Desktop/test.tre", required=True)
    args = parser.parse_args()

main(args.pathAndInfilename)

print colored("  --- End of Script ---", 'cyan')

