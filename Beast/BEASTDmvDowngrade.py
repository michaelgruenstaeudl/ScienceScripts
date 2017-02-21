#!/usr/bin/env python2
'''Downgrade output (dmv values) of BEAST'''
__author__ = "Michael Gruenstaeudl, PhD"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.07.10.1100"
__status__ = "Testing"

#########################
### IMPORT OPERATIONS ###
#########################

import numpy, sys
sys.path.insert(0,"/home/michael/git/ScienceScripts/")
import GeneralFileOperations as GFO
import GeneralStringOperations as GSO

###################
### DEFINITIONS ###
###################

def main():
    ind = GFO.loadRL(sys.argv[1])
    outd = []
    for line in ind:
        if "dmv={" not in line:
            outd.append(line)
        if "dmv={" in line:
            outd.append(line.split("[")[0])
            outd = outd + GSO.csplit(line,"[")[1:3]  # Do not outd.append()
            for listElem in GSO.csplit(line,"[")[3:]:
                cntr = 0
                for part in GSO.csplit(listElem,"]"):
                    if GSO.iseven(cntr):
                        dmvL = GSO.exstr(part,"{","}").split(",")
                        dmvL = [float(val) for val in dmvL]
                        outd.append("[&dmv="+str(numpy.mean(dmvL)))
                    else:
                        outd.append(part)
                    cntr =+ 1
            #outd.append(GSO.csplit(line,"[")[-1])
    #print(outd)
    outf = open(GSO.rmext(sys.argv[1])+".mod.trees","w")
    outf.write(''.join(outd)) # Do not '\n'.join
    outf.close()

###############
### EXECUTE ###
###############
main()

