#!/usr/bin/env python2
"""Beauti-Automator: Batch Generation of XML Input Files for BEAST and *BEAST."""
__author__ = "Michael Gruenstaeudl, PhD"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.06.12.0430"
__status__ = "Working"

#########################
### IMPORT OPERATIONS ###
#########################

import argparse, collections, re, sys
from Bio import SeqIO
from cStringIO import StringIO
from progress.bar import Bar
from termcolor import colored

########################
### GLOBAL VARIABLES ###
########################

glob_dict = {   "a": "pop1", "b": "pop2", "c": "pop3", "d": "pop4", "e": "pop5", "f": "pop6", 
                "g": "pop7", "h": "pop7", "i": "pop8", "o": "out4"  }
psd = "/home/michael/git/ScienceScripts/"                                       # Path to script directory

###################
### DEFINITIONS ###
###################

def generate_AR1(mode,counter,mydict):
    outlist = ['\t\t<sequence><taxon idref="'+key+'"/>'+str(mydict[key].seq)+'</sequence>' for key in mydict]
    if mode == "1":
        return '\t<alignment id="alignment" dataType="nucleotide">\n'+'\n'.join(outlist)+'\n\t</alignment>\n'
    if mode == "2":
        return '\t<alignment id="alignment'+str(counter)+'" dataType="nucleotide">\n'+'\n'.join(outlist)+'\n\t</alignment>\n'

def generate_OT1a(mydict):
    outlist = ['\t\t<taxon id="'+key+'"><attr name="species">'+glob_dict[key[0]]+'</attr></taxon>' for key in mydict]
    return '\n'.join(outlist)

def generate_OT1b(mydict):
    outlist = ['\t\t<taxon id="'+key+'"/>' for key in mydict]
    return '\n'.join(outlist)

def generate_OT2(mydict):
    unique_keystarts = set([key[0] for key in mydict])                          # set() gives unique values from list
    outstring = ""
    for letter in unique_keystarts:
        tmplist = ['\t\t\t<taxon idref="'+key+'"/>\n' for key in mydict if key[0]==letter]
        tmpstring = '\t\t<sp id="'+glob_dict[letter]+'">\n'+''.join(tmplist)+'\t\t</sp>\n'
        outstring += tmpstring
    return outstring

def generate_OT3(mydict):
    unique_keystarts = set([key[0] for key in mydict])
    outlist = ['\t\t<sp idref="'+glob_dict[letter]+'"/>' for letter in unique_keystarts]  # list comprehension
    return '\n'.join(outlist)

def generate_RE(mode, counter, pwd, fileprefix):
    handle = open(pwd+fileprefix+".txt").read()
    if mode == "1":
        return handle.replace("gene_NN.", "").replace("alignmentNN", "alignment")
    if mode == "2":
        return handle.replace("gene_NN", "gene"+str(counter)).replace("alignmentNN", "alignment"+str(counter))

def splitkeepsep(astring, sep):                                                 # splits a string by separator, but keeps separator; inspired by http://programmaticallyspeaking.com/
    return reduce(lambda acc, elem: acc + [elem] if elem == sep else acc[:-1] + [acc[-1] + elem], re.split("(%s)" % re.escape(sep), astring)[1:], [])

############
### MAIN ###
############

def main(mode, pwd, infilename, generations, options):

    logEvery = str(int(generations)/2000)                                       # setting up logEvery for *BEAST
    if logEvery < 1:
        logEvery = 1

    infile = open(pwd+infilename).read()

    if infile.count("#NEXUS") > 1:                                              # if infile contains multiple instances of "#NEXUS"
        alist = splitkeepsep(infile, "#NEXUS")                                  # split these instances, but keep the seperator (i.e. "#NEXUS")
        alist = filter(None,alist)                                              # remove all empty list elements
    else:
        alist = [infile]                                                        # else, create a list with just one element

    mylists = collections.defaultdict(list)                                     # initialization with a default dictionary

    if mode == "1":
        print "  Selected mode:",colored("BEAST","magenta")
        if options == "1":
            print "  Adding code for:",colored("Marginal Likelihoods","magenta")
        bar = Bar('  Generating XML files', max=len(alist))
        for nexuscounter,element in enumerate(alist,start=1):                   # loop through alist, also keep a counter
            handle = StringIO(element)                                          # convert string to file object
            mydict = SeqIO.to_dict(SeqIO.parse(handle, "nexus"))                # read file object and parse into dictionary
            filename = infilename[:-4]+".GeneTree.gene"+str(nexuscounter)

        # Loading Alignment Repeating (AR) elements; inside the loop
            mylists["AR1"] = generate_AR1(mode,nexuscounter, mydict)

        # Loading Repeating elements (RE); inside the loop
            for REnumber in range(1,24):
                mylists["RE"+str(REnumber)] = generate_RE(mode, nexuscounter, psd+"BeautiAutomatorRepository/", "RE"+str(REnumber))

        # Loading Stationary elements (SE)
            for SEnumber in range(1,16):
                mylists["SE"+str(SEnumber)] = open(psd+"BeautiAutomatorRepository/"+"SE"+str(SEnumber)+".txt").read()

            for SPnumber in range(1,4):
                mylists["SP"+str(SPnumber)] = open(psd+"BeautiAutomatorRepository/"+"SP"+str(SPnumber)+".txt").read()

        # Loading Marginal likelihood elements (MLE)
            mylists["MLE"] = open(psd+"BeautiAutomatorRepository/"+"MLE.txt").read()

            results = generateBEAST(mydict, mylists, filename, generations, logEvery, options)

            outfile = open(pwd+filename+".xml","w")
            outfile.write(results)                                              # outfile INSIDE the nexuscounter-loop
            outfile.close()
            bar.next()
        bar.finish()


    if mode == "2":
        print "  Selected mode:",colored("starBEAST","magenta")
        if options == "1":
            print "  Adding code for:",colored("Marginal Likelihoods","magenta")
        bar = Bar("  Generating XML files", max=len(alist))
        for nexuscounter,element in enumerate(alist,start=1):                   # loop through alist, also keep a counter
            handle = StringIO(element)                                          # convert string to file object
            mydict = SeqIO.to_dict(SeqIO.parse(handle, "nexus"))                # read file object and parse into dictionary
            
        # Loading Alignment Repeating (AR) elements; inside the loop
            mylists["AR1"].append(generate_AR1(mode,nexuscounter, mydict))

        # Loading Repeating elements (RE); inside the loop
            for REnumber in range(1,24):
                mylists["RE"+str(REnumber)].append(generate_RE(mode, nexuscounter, psd+"BeautiAutomatorRepository/", "RE"+str(REnumber)))

            bar.next()

        # Loading Stationary elements (SE); outside the loop
        for SEnumber in range(1,16):                                            # must not be in for-loop from above
            mylists["SE"+str(SEnumber)] = open(psd+"BeautiAutomatorRepository/"+"SE"+str(SEnumber)+".txt").read()

        # Loading Marginal likelihood elements (MLE); outside the loop
        mylists["MLE"] = open(psd+"BeautiAutomatorRepository/"+"MLE.txt").read()

        filename = infilename[:-4]+".SpeciesTree"                               # for starBeast, gene number is not part of the
        results = generateStarBEAST(mydict, mylists, filename, generations, logEvery, options)

        outfile = open(pwd+filename+".xml","w")                                 # outfile OUTSIDE the nexuscounter-loop
        outfile.write(results)
        outfile.close()
        bar.finish()


def generateBEAST(mydict, mylists, infilename, generations, logEvery, options):

    # Loading One-Time (OT) elements; outside the loop
    mylists["OT1b"] = generate_OT1b(mydict)

    # Replacements via RegEx
    mylists["SE3"] = re.sub(r'value="[^"]*"','value="0.0088"',mylists["SE3"])   # difference between BEAST and starBEAST

    start = mylists["SE10"].find('\t</operators>')
    end = mylists["SE10"].find('<prior id="prior">', start)+18
    mylists["SE10"] = mylists["SE10"][start:end]
    mylists["SE10"] = re.sub(r'chainLength="[^"]*"','chainLength="'+str(generations)+'"',mylists["SE10"])
    mylists["SE10"] = re.sub(r'operatorAnalysis="[^"]*.ops"','',mylists["SE10"])  

    start = mylists["SE11"].find('\t\t\t\t<oneOnXPrior>')
    end = mylists["SE11"].find('</oneOnXPrior>', start)+14
    mylists["SE11"] = mylists["SE11"][start:end]
    mylists["SE11"] = re.sub(r'idref="[^"]*"','idref="constant.popSize"',mylists["SE11"])

    end = mylists["SE12"].find('<column label="PopMean"')
    mylists["SE12"] = mylists["SE12"][:end]
    mylists["SE12"] = re.sub(r'logEvery="[^"]*"','logEvery="'+str(logEvery)+'"',mylists["SE12"])
    mylists["SE12"] = re.sub(r'<speciesCoalescent*','',mylists["SE12"])

    end = mylists["SE13"].find("<speciesCoalescent")
    mylists["SE13"] = mylists["SE13"][:end]
    mylists["SE13"] = re.sub(r'logEvery="[^"]*"','logEvery="'+str(logEvery)+'"',mylists["SE13"])
    mylists["SE13"] = re.sub(r'fileName="[^"]*.log"','fileName="'+infilename+'.log"',mylists["SE13"])

    mylists["RE23"] = re.sub(r'logEvery="[^"]*"','logEvery="'+str(logEvery)+'"',mylists["RE23"])
    mylists["RE23"] = re.sub(r'fileName="[^"]*.trees"','fileName="'+infilename+'.trees"',mylists["RE23"])

    mylists["MLE"] = re.sub(r'chainLength="[^"]*"','chainLength="'+str(generations)+'"',mylists["MLE"])
    mylists["MLE"] = re.sub(r'logEvery="[^"]*"','logEvery="'+str(logEvery)+'"',mylists["MLE"])
    mylists["MLE"] = re.sub(r'fileName="[^"]*.MargLikeEst.log"','fileName="'+infilename+'.MargLikeEst.log"',mylists["MLE"])

    results =   mylists["SE1"] +\
                mylists["OT1b"] +\
                mylists["SE2"] + '\n\n' +\
                mylists["AR1"] + '\n\n' +\
                mylists["RE1"] + '\n\n' +\
                mylists["SE3"] + '\n\n' +\
                mylists["RE2"] + '\n\n' +\
                mylists["RE3"] + '\n\n' +\
                mylists["SP1"] + '\n\n' +\
                mylists["RE4"] + '\n\n' +\
                mylists["RE5"] + '\n\n' +\
                mylists["RE6"] + '\n\n' +\
                mylists["SP2"] +\
                mylists["RE8"] + '\n\n' +\
                mylists["RE9"] + '\n\n' +\
                mylists["RE11"] + '\n\n' +\
                mylists["SP3"] +\
                mylists["RE12"] +\
                mylists["SE10"] + '\n\n' +\
                mylists["RE13"] + '\n\n' +\
                mylists["RE14"] + '\n\n' +\
                mylists["SE11"] + '\n\n' +\
                '\t\t\t\t<coalescentLikelihood idref="coalescent"/>\n\t\t\t</prior>\n\t\t\t<likelihood id="likelihood">\n' +\
                mylists["RE15"] +\
                mylists["SE12"] + '\n\n' +\
                mylists["RE16"] + '\n\n' +\
                mylists["RE17"] +\
                mylists["SE13"] +\
                mylists["RE18"] +\
                '\t\t\t\t<parameter idref="constant.popSize"/>\n' +\
                mylists["RE19"] +\
                mylists["RE20"] +\
                mylists["RE21"] +\
                mylists["RE22"] +\
                '\t\t\t<coalescentLikelihood idref="coalescent"/>\n\t\t\t</log>\n' +\
                mylists["RE23"] + '\n\n'                

    if options == "0":
        results = results + mylists["SE15"]  

    if options == "1":
        results = results + mylists["MLE"]

    return results


def generateStarBEAST(mydict, mylists, infilename, generations, logEvery, options):

    # Loading One-Time (OT) elements; outside the loop
    mylists["OT1a"] = generate_OT1a(mydict)
    mylists["OT2"] = generate_OT2(mydict)
    mylists["OT3"] = generate_OT3(mydict)

    # Replacements via RegEx
    mylists["SE10"] = re.sub(r'chainLength="[^"]*"','chainLength="'+str(generations)+'"',mylists["SE10"])
    mylists["SE10"] = re.sub(r'operatorAnalysis="[^"]*.ops"','operatorAnalysis="'+infilename+'.ops"',mylists["SE10"])

    mylists["SE12"] = re.sub(r'logEvery="[^"]*"','logEvery="'+str(logEvery)+'"',mylists["SE12"])

    mylists["SE13"] = re.sub(r'logEvery="[^"]*"','logEvery="'+str(logEvery)+'"',mylists["SE13"])
    mylists["SE13"] = re.sub(r'fileName="[^"]*.log"','fileName="'+infilename+'.log"',mylists["SE13"])

    mylists["SE14"] = re.sub(r'logEvery="[^"]*"','logEvery="'+str(logEvery)+'"',mylists["SE14"])
    mylists["SE14"] = re.sub(r'fileName="[^"]*.species.trees"','fileName="'+infilename+'.species.trees"',mylists["SE14"])

    mylists["RE23"] = re.sub(r'logEvery="[^"]*"','logEvery="'+str(logEvery)+'"','\n'.join(mylists["RE23"]))
    mylists["RE23"] = re.sub(r'fileName="[^"]*.gene','fileName="'+infilename+'.gene',mylists["RE23"])

    mylists["MLE"] = re.sub(r'chainLength="[^"]*"','chainLength="'+str(generations)+'"',mylists["MLE"])
    mylists["MLE"] = re.sub(r'logEvery="[^"]*"','logEvery="'+str(logEvery)+'"',mylists["MLE"])
    mylists["MLE"] = re.sub(r'fileName="[^"]*.MargLikeEst.log"','fileName="'+infilename+'.MargLikeEst.log"',mylists["MLE"])


    results =   mylists["SE1"] +\
                mylists["OT1a"] +\
                mylists["SE2"] + '\n\n' +\
                '\n'.join(mylists["AR1"]) + '\n\n' +\
                '\n'.join(mylists["RE1"]) + '\n\n' +\
                mylists["SE3"] + '\n\n' +\
                '\n'.join(mylists["RE2"]) + '\n\n' +\
                '\n'.join(mylists["RE3"]) + '\n\n' +\
                '\n'.join(mylists["RE4"]) + '\n\n' +\
                '\n'.join(mylists["RE5"]) + '\n\n' +\
                '\n'.join(mylists["RE6"]) + '\n\n' +\
                mylists["SE4"] +\
                mylists["OT2"] +\
                mylists["SE5"] +\
                '\n'.join(mylists["RE7"]) +\
                mylists["SE6"] +\
                mylists["OT3"] +\
                mylists["SE7"] + '\n\n' +\
                '\n'.join(mylists["RE8"]) + '\n\n' +\
                '\n'.join(mylists["RE9"]) + '\n\n' +\
                mylists["SE8"] + '\n\n' +\
                '\n'.join(mylists["RE10"]) + '\n\n' +\
                mylists["SE9"] + '\n\n' +\
                '\n'.join(mylists["RE11"]) + '\n\n' +\
                '\n'.join(mylists["RE12"]) + '\n\n' +\
                mylists["SE10"] + '\n\n' +\
                '\n'.join(mylists["RE13"]) + '\n\n' +\
                '\n'.join(mylists["RE14"]) + '\n\n' +\
                mylists["SE11"] + '\n\n' +\
                '\n'.join(mylists["RE15"]) + '\n\n' +\
                mylists["SE12"] + '\n\n' +\
                '\n'.join(mylists["RE16"]) + '\n\n' +\
                '\n'.join(mylists["RE17"]) +\
                mylists["SE13"] + '\n\n' +\
                '\n'.join(mylists["RE18"]) + '\n\n' +\
                '\n'.join(mylists["RE19"]) + '\n\n' +\
                '\n'.join(mylists["RE20"]) + '\n\n' +\
                '\n'.join(mylists["RE21"]) + '\n\n' +\
                '\n'.join(mylists["RE22"]) +\
                mylists["SE14"] + '\n\n' +\
                mylists["RE23"] + '\n\n'

    if options == "0":
        results = results + mylists["SE15"]  

    if options == "1":
        results = results + mylists["MLE"]

    return results

###############
### EXECUTE ###
###############

print ""
print colored("  Script name: "+sys.argv[0], 'cyan')
print colored("  Author: "+__author__, 'cyan')
print colored("  Version: "+__version__, 'cyan')
print ""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Batch Generation of XML Input Files for BEAST and *BEAST; 2014 Michael Gruenstaeudl')
    parser.add_argument('-m','--mode', help='1 = BEAST, 2 = starBEAST', default="1", required=True)
    parser.add_argument('-p','--pwd', help='/path/to/working/dir', default="~/Desktop/", required=True)
    parser.add_argument('-i','--infilename', help='name of input NEXUS file', default="test.nex", required=True)
    parser.add_argument('-g','--generations', help='number of MCMC generations', default="2000000", required=True)
    parser.add_argument('-o','--options', help='0 = none, 1 = Marginal Likelihood Estimation', default="0", required=False)
    args = parser.parse_args()

main(args.mode, args.pwd, args.infilename, args.generations, args.options)

print ""
print colored("  Done.", 'cyan')
print ""

