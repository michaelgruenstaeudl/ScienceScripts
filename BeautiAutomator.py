#!/usr/bin/env python2
'''BeautiAutomator: Batch Generation of XML Input Files for BEAST
and *BEAST.'''
__author__ = "Michael Gruenstaeudl, PhD"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.08.27.0300"
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


class makeOT1a(XMLtaglist):
    def __init__(self, inDict):
        super(makeOT1a, self).__init__()
        self.fromdict(inDict)

# Example output of TFL:
# <taxon id="b26"><attr name="species">pop2</attr></taxon>
    def build_node(self, key, val):
        el = ET.Element("taxon")               # Generates: <taxon >
        subel = ET.SubElement(el, "attr")      # Generates: <attr >
        subel.set("name", "species")           # Generates: name="species"
        subel.text = key[0].upper()            # Generates: B
        el.set("id", key)                      # Generates: id="b0026"
        return el

###############
# DEFINITIONS #
###############


def makeAR1(mode, cntr, inDict):
    outLst = ['\t\t<sequence><taxon idref="' + key + '"/>' +
              str(inDict[key]) + '</sequence>' for key in inDict]
    if mode == "1":
        return '\t<alignment id="alignment" dataType="nucleotide">\n' +\
               '\n'.join(outLst)+'\n\t</alignment>\n'
    if mode == "2":
        return '\t<alignment id="alignment' + cntr +\
               '" dataType="nucleotide">\n' + '\n'.join(outLst) +\
               '\n\t</alignment>\n'


def makeOT1b(inDict):
    outLst = ['\t\t<taxon id="' + key + '"/>' for key in inDict]
    return '\n'.join(outLst)


def makeOT2(inDict):
    # set() returns unique values from list
    uniqKeyStrt = set([key[0] for key in inDict])
    outStr = ""
    for ltr in uniqKeyStrt:
        # Using list comprehension in TFL
        tmpList = ['\t\t\t<taxon idref="' + key + '"/>\n'
                   for key in inDict if key[0] == ltr]
        tmpStr = '\t\t<sp id="' + ltr.upper() + '">\n' +\
                 ''.join(tmpList) + '\t\t</sp>\n'
        outStr += tmpStr
    return outStr


def makeOT3(inDict):
    uniqKeyStrt = set([key[0] for key in inDict])
    # Using list comprehension in TFL
    outLst = ['\t\t<sp idref="' + ltr.upper() + '"/>'
              for ltr in uniqKeyStrt]
    return '\n'.join(outLst)


def makeRE(mode, cntr, pwd, filePrefix):
    handle = open(pwd+filePrefix+".txt").read()
    if mode == "1":
        return (handle.replace("gene_NN.", "")
                .replace("alignmentNN", "alignment"))
    if mode == "2":
        return (handle.replace("gene_NN", "gene" + cntr)
                .replace("alignmentNN", "alignment" + cntr))


def makeModel(mode, modelDict, cntr, pwd, filePrefix):

    if not modelDict["nst"]:
        handle = open(pwd + filePrefix + ".txt").read()

    if modelDict["nst"]:

        if modelDict["base"] != "equal":
            # Calculate base frequencies
            baseFreqs = modelDict["base"].strip("(").rstrip(")").split()
            baseFreqs.append(str(1-sum([float(x) for x in baseFreqs])))
            baseFreqs = " ".join(baseFreqs)

        if modelDict["base"] == "equal":
            baseFreqs = "0.25 0.25 0.25 0.25"

        if modelDict["nst"] == "2":
            handle = open(pwd + filePrefix + "_HKY" + ".txt").read()
            # Replace base frequencies
            handle = handle.replace("baseFreqs", baseFreqs)
            # Replace Kappa values
            handle = handle.replace('kappa" value="2.0"',
                                    'kappa" value="' + modelDict["tratio"] + '"')
            # Adding site model info
            handle += '\t\t\t<HKYModel idref="gene_NN.hky"/>'

        if modelDict["nst"] == "6":
            handle = open(pwd + filePrefix + "_GTR" + ".txt").read()
            # Replace base frequencies
            handle = handle.replace("baseFreqs", baseFreqs)
            # Replace substitution rates
            substRates = modelDict["rmat"].strip("(").rstrip(")").split()
            for rate in substRates:
                handle = handle.replace('value="1.0"', 'value="'+rate+'"', 1)
            # Adding site model info
            handle += '\t\t\t<gtrModel idref="gene_NN.gtr"/>'

    if modelDict["pinv"]:
        handle += ''.join(['\n\t\t<proportionInvariant>\n',
                           '\t\t\t<parameter id="gene_NN.pInv" value="0.5" lower="0.0" upper="1.0"/>\n',
                           '\t\t</proportionInvariant>\n'])
        handle = handle.replace('pInv" value="0.5"',
                                'pInv" value="'+modelDict["pinv"]+'"')

    if modelDict["rates"] == "gamma":
        handle += ''.join(['\t\t<gammaShape gammaCategories="4">\n',
                           '\t\t\t<parameter id="gene_NN.alpha" value="0.5" lower="0.0"/>\n',
                           '\t\t</gammaShape>\n'])
        handle = handle.replace('alpha" value="0.5"',
                                'alpha" value="'+modelDict["shape"]+'"')

    # Closing site model info
    handle += '\n\t\t</substitutionModel>\n\t</siteModel>\n'

    if mode == "1":
        return (handle.replace("gene_NN.", "")
                .replace("alignmentNN", "alignment"))
    if mode == "2":
        return (handle.replace("gene_NN", "gene" + cntr)
                .replace("alignmentNN", "alignment" + cntr))


def generateBEAST(seqDict, myLists, inFn, nGens, logEvery):

    # Loading One-Time (OT) elements; outside the loop
    myLists["OT1b"] = makeOT1b(seqDict)

    # Replacements via RegEx
    # difference between BEAST and starBEAST
    myLists["SE3"] = re.sub(r'value="[^"]*"', 'value="0.0088"', myLists["SE3"])

    # SE10
    start = myLists["SE10"].find('\t</operators>')
    end = myLists["SE10"].find('<prior id="prior">', start) + 18
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
                             'fileName="' + inFn + '.log"',
                             myLists["SE13"])

    # RE23
    myLists["RE23"] = re.sub(r'logEvery="[^"]*"',
                             'logEvery="' + str(logEvery) + '"',
                             myLists["RE23"])

    myLists["RE23"] = re.sub(r'fileName="[^"]*.trees"',
                             'fileName="' + inFn + '.trees"',
                             myLists["RE23"])

#    # MLE
#    myLists["MLE"] = re.sub(r'chainLength="[^"]*"',
#                            'chainLength="' + str(nGens) + '"',
#                            myLists["MLE"])
#    myLists["MLE"] = re.sub(r'logEvery="[^"]*"',
#                            'logEvery="' + str(logEvery) + '"',
#                            myLists["MLE"])
#    myLists["MLE"] = re.sub(r'fileName="[^"]*.MargLikeEst.log"',
#                            'fileName="' + inFn + '.MargLikeEst.log"',
#                            myLists["MLE"])

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
        '\t\t\t\t<coalescentLikelihood idref="coalescent"/>\n' +\
        '\t\t\t</prior>\n' +\
        '\t\t\t<likelihood id="likelihood">\n' +\
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
        '\t\t\t<coalescentLikelihood idref="coalescent"/>\n' +\
        '\t\t\t</log>\n' +\
        myLists["RE23"] + '\n\n' +\
        myLists["SE15"]

#    if options == "0":
#        results = results + myLists["SE15"]
#    if options == "1":
#        results = results + myLists["MLE"]

    return results


def generateStarBEAST(seqDict, myLists, inFn, nGens, logEvery):

    # Only relevant for *BEAST (i.e. species tree)
    myLists["OT1a"] = makeOT1a(seqDict).get_str()

    myLists["OT2"] = makeOT2(seqDict)
    myLists["OT3"] = makeOT3(seqDict)

    # Replacements via RegEx

    # SE10
    myLists["SE10"] = re.sub(r'chainLength="[^"]*"',
                             'chainLength="' + str(nGens) + '"',
                             myLists["SE10"])
    myLists["SE10"] = re.sub(r'operatorAnalysis="[^"]*.ops"',
                             'operatorAnalysis="' + inFn + '.ops"',
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
                             'fileName="' + inFn + '.log"',
                             myLists["SE13"])

    # SE14
    myLists["SE14"] = re.sub(r'logEvery="[^"]*"',
                             'logEvery="' + str(logEvery) + '"',
                             myLists["SE14"])
    myLists["SE14"] = re.sub(r'fileName="[^"]*.species.trees"',
                             'fileName="' + inFn + '.species.trees"',
                             myLists["SE14"])

    # RE23
    myLists["RE23"] = re.sub(r'logEvery="[^"]*"',
                             'logEvery="' + str(logEvery) + '"',
                             '\n'.join(myLists["RE23"]))

    myLists["RE23"] = re.sub(r'fileName="[^"]*.gene',
                             'fileName="' + inFn + '.gene',
                             myLists["RE23"])

#    # MLE
#    myLists["MLE"] = re.sub(r'chainLength="[^"]*"',
#                            'chainLength="' + str(nGens) + '"',
#                            myLists["MLE"])
#    myLists["MLE"] = re.sub(r'logEvery="[^"]*"',
#                            'logEvery="' + str(logEvery) + '"',
#                            myLists["MLE"])
#    myLists["MLE"] = re.sub(r'fileName="[^"]*.MargLikeEst.log"',
#                            'fileName="' + inFn + '.MargLikeEst.log"',
#                            myLists["MLE"])

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
        myLists["RE23"] + '\n\n' +\
        myLists["SE15"]

#    if options == "0":
#        results = results + myLists["SE15"]
#    if options == "1":
#        results = results + myLists["MLE"]

    return results


def addLeadZeros(inDict):
    outDict = {}
    for key in inDict:
        newKey = key[0] + str(int(key[1:])).zfill(4)
        outDict[newKey] = inDict[key].seq
    return outDict


def extractModelInfo(inStr):
    # 1. Extract lines that starts with "lset"
    found = re.search("lset(.+?)\n", inStr).group(1)
    if found:
        modelDict = {}
        # EXAMPLE: "set crit=like; lset clock=no nst=6 base=(0.21901 0.26151 0.36629) rmat=(1.000000 5.187448 0.293936 0.293936 5.187448) rclass=(a b c c b a) rates=equal pinv=0.832102;"
        # remove trailing ";"
        tmpStr = found.rstrip(";")
        # remove info on "rclass"
        tmpStr = re.sub(r'rclass=(.+?)\)', '', tmpStr)
        # http://stackoverflow.com/questions/25518647/re-split-unless-following-character-an-integer-in-python/25518698#25518698
        aList = filter(None, re.split(r' (?!\d)', tmpStr))
        for elem in aList:
            bList = re.split('=', elem)
            modelDict[bList[0]] = bList[1]
    return modelDict


########
# MAIN #
########


def main(mode, pwd, inFn, nGens):

    # setting up logEvery for BEAST/*BEAST
    logEvery = str(int(nGens)/5000)

    infile = open(pwd + "/" + inFn).read()
    # if infile contains multiple instances of "#NEXUS"
    if infile.count("#NEXUS") > 1:
        # split these instances, but keep the seperator (i.e. "#NEXUS")
        nxsList = GSO.splitkeepsep(infile, "#NEXUS")
        # remove all empty list elements
        nxsList = filter(None, nxsList)
    else:
        # else, create a list with just one element
        nxsList = [infile]

    # Set up inFn without extension
    inFnStem = GSO.rmext(inFn)

    # initialization with a default dictionary
    myLists = collections.defaultdict(list)

    if mode == "1":
        print "\tSelected mode: " + colored("BEAST", "magenta")
#        if options == "1":
#            print "\tAdding code for:" +\
#                colored("Marginal Likelihoods", "magenta")

        # Showing progress bar
        bar = Bar('\tGenerating XML files', max=len(nxsList))

        # Loop through nxsList, also keep a cntr
        for nxsCntr, elem in enumerate(nxsList, start=1):

            # Extract model information
            modelDict = extractModelInfo(elem)

            # Convert string to file object
            handle = StringIO(elem)
            # Read file object and parse into dictionary
            inDict = SeqIO.to_dict(SeqIO.parse(handle, "nexus"))

            # Convert keys in inDict to appropriate numbering
            seqDict = addLeadZeros(inDict)

            # Loading Alignment Repeating (AR) elements; inside the loop
            myLists["AR1"] = makeAR1(mode, str(nxsCntr).zfill(3), seqDict)

            # Loading Repeating elements (RE); inside the loop
            for REnmbr in range(1, 5)+range(6, 24):
                myLists["RE" + str(REnmbr)] = makeRE(mode,
                                                     str(nxsCntr).zfill(3),
                                                     psd + "BARepo/",
                                                     "RE" + str(REnmbr))

            # Special RE case: Model info
            myLists["RE5"] = makeModel(mode, modelDict, str(nxsCntr).zfill(3),
                                       psd + "BARepo/", "RE5")

            # Loading Stationary elements (SE)
            for SEnmbr in range(1, 16):
                myLists["SE" + str(SEnmbr)] = open(psd + "BARepo/" +
                                                   "SE" + str(SEnmbr) +
                                                   ".txt").read()

            for SPnmbr in range(1, 4):
                myLists["SP" + str(SPnmbr)] = open(psd + "BARepo/" +
                                                   "SP" + str(SPnmbr) +
                                                   ".txt").read()

#            # Loading Marginal likelihood elements (MLE)
#            myLists["MLE"] = open(psd+"BARepo/"+"MLE.txt").read()

            fName = inFnStem + ".GeneTree.gene" + str(nxsCntr).zfill(3)

            results = generateBEAST(seqDict, myLists,
                                    fName, nGens, logEvery)

            # outF INSIDE the nxsCntr-loop
            outF = open(pwd + "/" + fName + ".xml", "w")
            outF.write(results)
            outF.close()
            bar.next()
        bar.finish()

    if mode == "2":
        print "\tSelected mode: " + colored("starBEAST", "magenta")
#        if options == "1":
#            print "    Adding code for:" +\
#                colored("Marginal Likelihoods", "magenta")

        bar = Bar("\tGenerating XML files", max=len(nxsList))
        # loop through nxsList, also keep a cntr
        for nxsCntr, elem in enumerate(nxsList, start=1):

            # Extract model information
            modelDict = extractModelInfo(elem)

            # convert string to file object
            handle = StringIO(elem)
            # read file object and parse into dictionary
            inDict = SeqIO.to_dict(SeqIO.parse(handle, "nexus"))

            # Convert keys in inDict to appropriate numbering
            seqDict = addLeadZeros(inDict)

            # Loading Alignment Repeating (AR) elements; inside the loop
            myLists["AR1"].append(makeAR1(mode,
                                          str(nxsCntr).zfill(3),
                                          seqDict))

            # Loading Repeating elements (RE); inside the loop
            for REnmbr in range(1, 5)+range(6, 24):
                myLists["RE" + str(REnmbr)].append(
                    makeRE(mode,
                           str(nxsCntr).zfill(3),
                           psd + "BARepo/",
                           "RE" + str(REnmbr)))

            # Special RE case: Model info
            myLists["RE5"].append(
                makeModel(mode,
                          modelDict,
                          str(nxsCntr).zfill(3),
                          psd + "BARepo/", "RE5"))

            bar.next()

        # Loading Stationary elements (SE); outside the loop
        # must not be in for-loop from above
        for SEnmbr in range(1, 16):
            myLists["SE" + str(SEnmbr)] = open(psd + "BARepo/" +
                                               "SE" + str(SEnmbr) +
                                               ".txt").read()

#        # Loading Marginal likelihood elements (MLE); outside the loop
#        myLists["MLE"] = open(psd + "BARepo/" + "MLE.txt").read()

        fName = inFnStem + ".SpeciesTree"

        results = generateStarBEAST(seqDict, myLists,
                                    fName, nGens, logEvery)

        # outF OUTSIDE the nxsCntr-loop
        outF = open(pwd + "/" + fName + ".xml", "w")
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
                        '--infile',
                        help='name of input NEXUS file',
                        default="test.nex",
                        required=True)

    parser.add_argument('-g',
                        '--gens',
                        help='number of MCMC generations',
                        default="50000000",
                        required=True)

#    parser.add_argument('-o',
#                        '--options',
#                        help='0 = none, 1 = Marginal Likelihood Estimation',
#                        default="0",
#                        required=False)

    args = parser.parse_args()

main(args.mode, args.pwd, args.infile, args.gens)

print ""
print colored("    Done.", 'cyan')
print ""
