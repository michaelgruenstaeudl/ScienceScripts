#!/usr/bin/env python2
'''Gene tree simulation under the coalescent using DendroPy'''
__author__ = "Michael Gruenstaeudl, PhD"
__email__ = "gruenstaeudl.1@osu.edu"
__version__ = "2014.01.05.1200"
__status__ = "Working"


#####################
# IMPORT OPERATIONS #
#####################

import argparse
import dendropy
import numpy
import sys
from termcolor import colored


########
# MAIN #
########

def main(n_sp, n_loci, alleles_per_sp, max_age, output_type):

## STEP 1: Simulate a population/species tree
## STEP 1a: Set up a taxon set
    d = dendropy.DataSet()
    d.attach_taxon_set()
    tn_list = []
    for s in range(0, n_sp):
        letter = chr(s + ord("a"))
        tn_str = ">" + letter + "\nA"
        tn_list.append(tn_str)

    taxon_names = "\n".join(tn_list)
    d.read_from_string(taxon_names, "dnafasta")

## STEP 1B: Set up age vector
    ages = list(numpy.random.uniform(0, max_age, size=n_sp-1))
    # Sort the ages (this is critical to get an independent sample!)
    ages = sorted(ages, key=int)
    ages.append(max_age)

## STEP 1c: Set up population sizes
    ps = [1] * (n_sp * 2 + 1)

## STEP 1d: Simulate the species tree
    sp_tree = dendropy.treesim.pop_gen_tree(taxon_set=d.taxon_sets[0], ages=ages, pop_sizes=ps)
    if output_type == 2|3:
        print("# SPECIES TREE")
        print(sp_tree.as_newick_string() + ";")

## STEP 2: Simulate a set of gene trees
## STEP 2a: Set up the TaxonSetMapping object
    tsm_object = dendropy.TaxonSetMapping.create_contained_taxon_mapping(containing_taxon_set=d.taxon_sets[0], num_contained=alleles_per_sp)

## STEP 2b: Simulate the gene trees
    g_trees = []
    for t in range(0, n_loci):
        g_tree = dendropy.treesim.contained_coalescent(sp_tree, tsm_object)
        g_trees.append(g_tree.as_newick_string())

## STEP 3: Save gene trees to file
#    outFn = "geneTrees_SimUnderCoalescent_%sspecies_%sloci_%sallelesPerSp" % (n_sp, n_loci, alleles_per_sp)
#    outfile = open(outFn + ".nwk", "w")
#    out_str = ";\n".join(g_trees) + ";"
#    out_str = out_str.replace("_", "")
#    outfile.write(out_str)
#    outfile.close()

## STEP 3: Print gene trees to screen
    if output_type == 1|3:
        print("# GENE TREES")
        gt_out_str = ";\n".join(g_trees) + ";"
        gt_out_str = gt_out_str.replace("_", "")
        print(gt_out_str)

###########
# EXECUTE #
###########

#print ""
#print colored("    Script name: " + sys.argv[0], 'cyan')
#print colored("    Author: " + __author__, 'cyan')
#print colored("    Version: " + __version__, 'cyan')
#print ""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Gene tree simulation under the coalescent \
                     using DendroPy; 2015 Michael Gruenstaeudl; \
                     Example Usage: python2 ~/git/ScienceScripts/gtCoalSim.py > ~/Desktop/geneTrees.tre')

    parser.add_argument('-s',
                        '--n_sp',
                        help='number of species',
                        default=6,
                        required=False)

    parser.add_argument('-l',
                        '--n_loci',
                        help='number of loci',
                        default=5,
                        required=False)

    parser.add_argument('-p',
                        '--alleles_per_sp',
                        help='number of alleles per species',
                        default=10,
                        required=False)

    parser.add_argument('-a',
                        '--max_age',
                        help='age of the oldest divergence (in generations)',
                        default=40,
                        required=False)

    parser.add_argument('-o',
                        '--output_type',
                        help='1 = gene trees only; 2 = species tree only; 3 = gene and species trees',
                        default=3,
                        required=False)

    args = parser.parse_args()

main(int(args.n_sp), int(args.n_loci), int(args.alleles_per_sp), int(args.max_age), int(args.output_type))

#print ""
#print colored("    Done.", 'cyan')
#print ""
