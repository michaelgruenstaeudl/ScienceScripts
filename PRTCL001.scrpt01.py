import dendropy
import GeneralStringOperations as GSO
import sys

print ""
print "USAGE: <this_script.py> <msCmd> <name_of_treefile>"
print ""

## Get pop info from commandline
msCmd = sys.argv[2]
# extract relevant section
popLst = GSO.exstr(msCmd, "-I", "-")
# remove leading and trailing whitespaces, make into list
popLst = popLst.lstrip().rstrip().split()
# remove:   first element (merely number of pops, not size)
#           last element: the outgroup
popLst = popLst[1:-1]

## Generate dictionary for label replacement
lttrs = map(chr, range(97, 123))
# remove "o" from list
lttrs.remove("o")
# Initializations
Sum = 0
aDict = {}
# loop over popLst
for cntr, pop in enumerate(popLst, start=1):
    for num in range(Sum+1, Sum+int(pop)+1):
        aDict[num] = lttrs[cntr] + str(num).zfill(4)
    Sum += int(pop)
# include an outgroup
    aDict[Sum+1] = "o" + str(Sum+1).zfill(4)

print aDict

### Load trees and perform label replacement
#trees = dendropy.TreeList.get_from_path('outfile.tre', 'newick')
#for tree in trees:
#    for leaf in tree.leaf_nodes():
#        leaf.taxon.label = aDict[leaf.taxon.label]

### Save trees back to file
