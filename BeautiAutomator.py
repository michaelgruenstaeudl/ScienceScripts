#!/usr/bin/env python2
'''BeautiAutomator: Batch Generation of XML Input Files for BEAST
and *BEAST.'''
__author__ = "Michael Gruenstaeudl, PhD"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.08.25.1200"
__status__ = "Working"

#####################
# IMPORT OPERATIONS #
#####################

import argparse
import collections
import GeneralStringOperations as GSO
import re
import sys
import xml.etree.ElementTree as ET
from Bio import SeqIO
from cStringIO import StringIO
from progress.bar import Bar
from termcolor import colored

####################
# GLOBAL VARIABLES #
####################

# Path to script directory
psd = "/home/michael/git/ScienceScripts/"

###########
# CLASSES #
###########

# CURRENTLY UNUSED:
#class XMLgenerator:
#    def __init__(self):
#        self.tree = ET.ElementTree()
#
#    def get_str(self):
#        return self.tree.to_str()
#
#    def get_tree(self):
#        return self.tree


class XMLtaglist(object):
    def __init__(self):
        self.tags = []

    def fromdict(self, valOf):
        self.tags = [self.build_node(key, valOf[key]) for key in valOf]

    def get_str(self):
        return "\n".join([ET.tostring(el) for el in self.tags])

    def get_list(self):
        return self.tags


class OT1aGenerator(XMLtaglist):
    def __init__(self, myDict):
        super(OT1aGenerator, self).__init__()
        self.fromdict(myDict)

# Example output of TFL:
# <taxon id="b26"><attr name="species">pop2</attr></taxon>
    def build_node(self, key, val):
        el = ET.Element("taxon")               # Generates: <taxon >
        subel = ET.SubElement(el, "attr")      # Generates: <attr >
        subel.set("name", "species")           # Generates: name="species"
        subel.text = key[0].upper()            # Generates: B
        el.set("id", key)                      # Generates: id="b26"
        return el

###############
# DEFINITIONS #
###############


def mke_AR1(mode, cntr, myDict):
    outLst = ['\t\t<sequence><taxon idref="'+key+'"/>' +
              str(myDict[key].seq) + '</sequence>' for key in myDict]
    if mode == "1":
        return '\t<alignment id="alignment" dataType="nucleotide">\n' +\
               '\n'.join(outLst)+'\n\t</alignment>\n'
    if mode == "2":
        return '\t<alignment id="alignment' + str(cntr) +\
               '" dataType="nucleotide">\n' + '\n'.join(outLst) +\
               '\n\t</alignment>\n'

# LEGACYCODE:
#def mke_OT1a(myDict):
#    outLst = ['\t\t<taxon id="'+key+'"><attr name="species">'+glob_dict[key[0]]+'</attr></taxon>' for key in myDict]
#    return '\n'.join(outLst)


def mke_OT1b(myDict):
    outLst = ['\t\t<taxon id="'+key+'"/>' for key in myDict]
    return '\n'.join(outLst)


def mke_OT2(myDict):
    # set() returns unique values from list
    unique_keystarts = set([key[0] for key in myDict])
    outStr = ""
    for ltr in unique_keystarts:
        # Using list comprehension in TFL
        tmpList = ['\t\t\t<taxon idref="' + key + '"/>\n'
                   for key in myDict if key[0] == ltr]
        tmpStr = '\t\t<sp id="' + ltr.upper() + '">\n' +\
                 ''.join(tmpList) + '\t\t</sp>\n'
        outStr += tmpStr
    return outStr


def mke_OT3(myDict):
    unique_keystarts = set([key[0] for key in myDict])
    # Using list comprehension in TFL
    outLst = ['\t\t<sp idref="' + ltr.upper() + '"/>'
              for ltr in unique_keystarts]
    return '\n'.join(outLst)


def mke_RE(mode, cntr, pwd, filePrefix):
    handle = open(pwd+filePrefix+".txt").read()
    if mode == "1":
        return (handle.replace("gene_NN.", "")
                .replace("alignmentNN", "alignment"))
    if mode == "2":
        return (handle.replace("gene_NN", "gene"+str(cntr))
                .replace("alignmentNN", "alignment"+str(cntr)))


def generateBEAST(myDict, myLists, inFname, nGens, logEvery, options):

    # Loading One-Time (OT) elements; outside the loop
    myLists["OT1b"] = mke_OT1b(myDict)

    # Replacements via RegEx
    # difference between BEAST and starBEAST
    myLists["SE3"] = re.sub(r'value="[^"]*"', 'value="0.0088"', myLists["SE3"])

    # SE10
    start = myLists["SE10"].find('\t</operators>')
    end = myLists["SE10"].find('<prior id="prior">', start)+18
    myLists["SE10"] = myLists["SE10"][start:end]
    myLists["SE10"] = re.sub(r'chainLength="[^"]*"',
                             'chainLength="' + str(nGens) + '"',
                             myLists["SE10"])
    myLists["SE10"] = re.sub(r'operatorAnalysis="[^"]*.ops"',
                             '',
                             myLists["SE10"])

    # SE11
    start = myLists["SE11"].find('\t\t\t\t<oneOnXPrior>')
    end = myLists["SE11"].find('</oneOnXPrior>', start) + 14
    myLists["SE11"] = myLists["SE11"][start:end]
    myLists["SE11"] = re.sub(r'idref="[^"]*"',
                             'idref="constant.popSize"',
                             myLists["SE11"])

    # SE12
    end = myLists["SE12"].find('<column label="PopMean"')
    myLists["SE12"] = myLists["SE12"][:end]
    myLists["SE12"] = re.sub(r'logEvery="[^"]*"',
                             'logEvery="' + str(logEvery) + '"',
                             myLists["SE12"])
    myLists["SE12"] = re.sub(r'<speciesCoalescent*',
                             '',
                             myLists["SE12"])

    # SE13
    end = myLists["SE13"].find("<speciesCoalescent")
    myLists["SE13"] = myLists["SE13"][:end]
    myLists["SE13"] = re.sub(r'logEvery="[^"]*"',
                             'logEvery="' + str(logEvery) + '"',
                             myLists["SE13"])
    myLists["SE13"] = re.sub(r'fileName="[^"]*.log"',
                             'fileName="' + inFname + '.log"',
                             myLists["SE13"])

    # RE23
    myLists["RE23"] = re.sub(r'logEvery="[^"]*"',
                             'logEvery="' + str(logEvery) + '"',
                             myLists["RE23"])
    myLists["RE23"] = re.sub(r'fileName="[^"]*.trees"',
                             'fileName="'+inFname+'.trees"',
                             myLists["RE23"])

    # MLE
    myLists["MLE"] = re.sub(r'chainLength="[^"]*"',
                            'chainLength="' + str(nGens) + '"',
                            myLists["MLE"])
    myLists["MLE"] = re.sub(r'logEvery="[^"]*"',
                            'logEvery="' + str(logEvery) + '"',
                            myLists["MLE"])
    myLists["MLE"] = re.sub(r'fileName="[^"]*.MargLikeEst.log"',
                            'fileName="' + inFname + '.MargLikeEst.log"',
                            myLists["MLE"])

    results = myLists["SE1"] +\
        myLists["OT1b"] +\
        myLists["SE2"] + '\n\n' +\
        myLists["AR1"] + '\n\n' +\
        myLists["RE1"] + '\n\n' +\
        myLists["SE3"] + '\n\n' +\
        myLists["RE2"] + '\n\n' +\
        myLists["RE3"] + '\n\n' +\
        myLists["SP1"] + '\n\n' +\
        myLists["RE4"] + '\n\n' +\
        myLists["RE5"] + '\n\n' +\
        myLists["RE6"] + '\n\n' +\
        myLists["SP2"] +\
        myLists["RE8"] + '\n\n' +\
        myLists["RE9"] + '\n\n' +\
        myLists["RE11"] + '\n\n' +\
        myLists["SP3"] +\
        myLists["RE12"] +\
        myLists["SE10"] + '\n\n' +\
        myLists["RE13"] + '\n\n' +\
        myLists["RE14"] + '\n\n' +\
        myLists["SE11"] + '\n\n' +\
        '\t\t\t\t<coalescentLikelihood idref="coalescent"/>' +\
        '\n\t\t\t</prior>\n\t\t\t<likelihood id="likelihood">\n' +\
        myLists["RE15"] +\
        myLists["SE12"] + '\n\n' +\
        myLists["RE16"] + '\n\n' +\
        myLists["RE17"] +\
        myLists["SE13"] +\
        myLists["RE18"] +\
        '\t\t\t\t<parameter idref="constant.popSize"/>\n' +\
        myLists["RE19"] +\
        myLists["RE20"] +\
        myLists["RE21"] +\
        myLists["RE22"] +\
        '\t\t\t<coalescentLikelihood idref="coalescent"/>' +\
        '\n\t\t\t</log>\n' +\
        myLists["RE23"] + '\n\n'

    if options == "0":
        results = results + myLists["SE15"]

    if options == "1":
        results = results + myLists["MLE"]

    return results


def generateStarBEAST(myDict, myLists, inFname, nGens, logEvery, options):

    # Only relevant for *BEAST (i.e. species tree)
    myLists["OT1a"] = OT1aGenerator(myDict).get_str()

    myLists["OT2"] = mke_OT2(myDict)
    myLists["OT3"] = mke_OT3(myDict)

    # Replacements via RegEx

    # SE10
    myLists["SE10"] = re.sub(r'chainLength="[^"]*"',
                             'chainLength="' + str(nGens) + '"',
                             myLists["SE10"])
    myLists["SE10"] = re.sub(r'operatorAnalysis="[^"]*.ops"',
                             'operatorAnalysis="' + inFname + '.ops"',
                             myLists["SE10"])

    # SE12
    myLists["SE12"] = re.sub(r'logEvery="[^"]*"',
                             'logEvery="' + str(logEvery) + '"',
                             myLists["SE12"])

    # SE13
    myLists["SE13"] = re.sub(r'logEvery="[^"]*"',
                             'logEvery="' + str(logEvery) + '"',
                             myLists["SE13"])
    myLists["SE13"] = re.sub(r'fileName="[^"]*.log"',
                             'fileName="' + inFname + '.log"',
                             myLists["SE13"])

    # SE14
    myLists["SE14"] = re.sub(r'logEvery="[^"]*"',
                             'logEvery="' + str(logEvery) + '"',
                             myLists["SE14"])
    myLists["SE14"] = re.sub(r'fileName="[^"]*.species.trees"',
                             'fileName="' + inFname + '.species.trees"',
                             myLists["SE14"])

    # RE23
    myLists["RE23"] = re.sub(r'logEvery="[^"]*"',
                             'logEvery="' + str(logEvery) + '"',
                             '\n'.join(myLists["RE23"]))
    myLists["RE23"] = re.sub(r'fileName="[^"]*.gene',
                             'fileName="' + inFname + '.gene',
                             myLists["RE23"])

    # MLE
    myLists["MLE"] = re.sub(r'chainLength="[^"]*"',
                            'chainLength="' + str(nGens) + '"',
                            myLists["MLE"])
    myLists["MLE"] = re.sub(r'logEvery="[^"]*"',
                            'logEvery="' + str(logEvery) + '"',
                            myLists["MLE"])
    myLists["MLE"] = re.sub(r'fileName="[^"]*.MargLikeEst.log"',
                            'fileName="' + inFname + '.MargLikeEst.log"',
                            myLists["MLE"])

    results = myLists["SE1"] +\
        myLists["OT1a"] +\
        myLists["SE2"] + '\n\n' +\
        '\n'.join(myLists["AR1"]) + '\n\n' +\
        '\n'.join(myLists["RE1"]) + '\n\n' +\
        myLists["SE3"] + '\n\n' +\
        '\n'.join(myLists["RE2"]) + '\n\n' +\
        '\n'.join(myLists["RE3"]) + '\n\n' +\
        '\n'.join(myLists["RE4"]) + '\n\n' +\
        '\n'.join(myLists["RE5"]) + '\n\n' +\
        '\n'.join(myLists["RE6"]) + '\n\n' +\
        myLists["SE4"] +\
        myLists["OT2"] +\
        myLists["SE5"] +\
        '\n'.join(myLists["RE7"]) +\
        myLists["SE6"] +\
        myLists["OT3"] +\
        myLists["SE7"] + '\n\n' +\
        '\n'.join(myLists["RE8"]) + '\n\n' +\
        '\n'.join(myLists["RE9"]) + '\n\n' +\
        myLists["SE8"] + '\n\n' +\
        '\n'.join(myLists["RE10"]) + '\n\n' +\
        myLists["SE9"] + '\n\n' +\
        '\n'.join(myLists["RE11"]) + '\n\n' +\
        '\n'.join(myLists["RE12"]) + '\n\n' +\
        myLists["SE10"] + '\n\n' +\
        '\n'.join(myLists["RE13"]) + '\n\n' +\
        '\n'.join(myLists["RE14"]) + '\n\n' +\
        myLists["SE11"] + '\n\n' +\
        '\n'.join(myLists["RE15"]) + '\n\n' +\
        myLists["SE12"] + '\n\n' +\
        '\n'.join(myLists["RE16"]) + '\n\n' +\
        '\n'.join(myLists["RE17"]) +\
        myLists["SE13"] + '\n\n' +\
        '\n'.join(myLists["RE18"]) + '\n\n' +\
        '\n'.join(myLists["RE19"]) + '\n\n' +\
        '\n'.join(myLists["RE20"]) + '\n\n' +\
        '\n'.join(myLists["RE21"]) + '\n\n' +\
        '\n'.join(myLists["RE22"]) +\
        myLists["SE14"] + '\n\n' +\
        myLists["RE23"] + '\n\n'

    if options == "0":
        results = results + myLists["SE15"]
    if options == "1":
        results = results + myLists["MLE"]

    return results


########
# MAIN #
########


def main(mode, pwd, inFname, nGens, options):

    # setting up logEvery for BEAST/*BEAST
    logEvery = str(int(nGens)/2000)

    infile = open(pwd+inFname).read()

    # if infile contains multiple instances of "#NEXUS"
    if infile.count("#NEXUS") > 1:
        # split these instances, but keep the seperator (i.e. "#NEXUS")
        alist = GSO.splitkeepsep(infile, "#NEXUS")
        # remove all empty list elements
        alist = filter(None, alist)
    else:
        # else, create a list with just one element
        alist = [infile]

    # initialization with a default dictionary
    myLists = collections.defaultdict(list)

    if mode == "1":
        print "    Selected mode:" + colored("BEAST", "magenta")
        if options == "1":
            print "    Adding code for:" +\
                colored("Marginal Likelihoods", "magenta")
        # Showing progress bar
        bar = Bar('    Generating XML files', max=len(alist))
        # Loop through alist, also keep a cntr
        for nxsCntr, element in enumerate(alist, start=1):
            # convert string to file object
            handle = StringIO(element)
            # read file object and parse into dictionary
            myDict = SeqIO.to_dict(SeqIO.parse(handle, "nexus"))
            filename = inFname[:-4]+".GeneTree.gene"+str(nxsCntr)

        # Loading Alignment Repeating (AR) elements; inside the loop
            myLists["AR1"] = mke_AR1(mode, nxsCntr, myDict)

        # Loading Repeating elements (RE); inside the loop
            for REnmbr in range(1, 24):
                myLists["RE" + str(REnmbr)] = mke_RE(mode, nxsCntr,
                                                     psd + "BARepo/",
                                                     "RE" + str(REnmbr))

        # Loading Stationary elements (SE)
            for SEnmbr in range(1, 16):
                myLists["SE" + str(SEnmbr)] = open(psd + "BARepo/" +
                                                   "SE" + str(SEnmbr) +
                                                   ".txt").read()

            for SPnmbr in range(1, 4):
                myLists["SP" + str(SPnmbr)] = open(psd + "BARepo/" +
                                                   "SP" + str(SPnmbr) +
                                                   ".txt").read()

        # Loading Marginal likelihood elements (MLE)
            myLists["MLE"] = open(psd+"BARepo/"+"MLE.txt").read()

            results = generateBEAST(myDict, myLists, filename,
                                    nGens, logEvery, options)

            # outF INSIDE the nxsCntr-loop
            outF = open(pwd+filename+".xml", "w")
            outF.write(results)
            outF.close()
            bar.next()
        bar.finish()

    if mode == "2":
        print "    Selected mode:", colored("starBEAST", "magenta")
        if options == "1":
            print "    Adding code for:" +\
                colored("Marginal Likelihoods", "magenta")
        bar = Bar("    Generating XML files", max=len(alist))
        # loop through alist, also keep a cntr
        for nxsCntr, element in enumerate(alist, start=1):
            # convert string to file object
            handle = StringIO(element)
            # read file object and parse into dictionary
            myDict = SeqIO.to_dict(SeqIO.parse(handle, "nexus"))

            # Loading Alignment Repeating (AR) elements; inside the loop
            myLists["AR1"].append(mke_AR1(mode, nxsCntr, myDict))

            # Loading Repeating elements (RE); inside the loop
            for REnmbr in range(1, 24):
                myLists["RE"+str(REnmbr)].append(mke_RE(mode, nxsCntr,
                                                        psd + "BARepo/",
                                                        "RE" + str(REnmbr)))

            bar.next()

        # Loading Stationary elements (SE); outside the loop
        # must not be in for-loop from above
        for SEnmbr in range(1, 16):
            myLists["SE"+str(SEnmbr)] = open(psd + "BARepo/" +
                                             "SE" + str(SEnmbr) +
                                             ".txt").read()

        # Loading Marginal likelihood elements (MLE); outside the loop
        myLists["MLE"] = open(psd + "BARepo/" + "MLE.txt").read()

        filename = GSO.rmext(inFname) + ".SpeciesTree"
        results = generateStarBEAST(myDict, myLists, filename,
                                    nGens, logEvery, options)

        # outF OUTSIDE the nxsCntr-loop
        outF = open(pwd+filename+".xml", "w")
        outF.write(results)
        outF.close()
        bar.finish()

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
        description='Batch Generation of XML Input Files for \
                     BEAST and *BEAST; 2014 Michael Gruenstaeudl')

    parser.add_argument('-m',
                        '--mode',
                        help='1 = BEAST, 2 = starBEAST',
                        default="1",
                        required=True)

    parser.add_argument('-p',
                        '--pwd',
                        help='/path/to/working/dir',
                        default="~/Desktop/",
                        required=True)

    parser.add_argument('-i',
                        '--infilename',
                        help='name of input NEXUS file',
                        default="test.nex",
                        required=True)

    parser.add_argument('-g',
                        '--generations',
                        help='number of MCMC generations',
                        default="50000000",
                        required=True)

    parser.add_argument('-o',
                        '--options',
                        help='0 = none, 1 = Marginal Likelihood Estimation',
                        default="0",
                        required=False)

    args = parser.parse_args()

main(args.mode, args.pwd, args.infilename, args.generations, args.options)

print ""
print colored("    Done.", 'cyan')
print ""
