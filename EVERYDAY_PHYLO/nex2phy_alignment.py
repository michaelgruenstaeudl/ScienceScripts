#!/usr/bin/env python2.7
import sys
from Bio import AlignIO
inFn = sys.argv[1]
inp = open(inFn, 'rU')
outp = open(inFn+'.phy', 'w')
aln = AlignIO.parse(inp, 'nexus')
AlignIO.write(aln, outp, 'phylip')
outp.close()
inp.close()

