#!/bin/bash

# Generate ReplaceDict from .NEX-file
INF=Nymph-IJPS-trnTF-inhot-mar_plus-genomes_11-2-17_exhot.nex
OUTF=${INF%.nex*}__ReplaceDict.txt
cat $INF | sed '/MATRIX/,$!d' | sed '/END\;/q' | awk '{print $1}' > $OUTF

# Conduct replacement
INF=Nymph-IJPS-trnTF-inhot-mar_plus-genomes_11-2-17_exhot.nex
IND=TreePosition_Nmexicana__ReplaceDict.txt
python2 /home/michael_science/git/michaelgruenstaeudl_ScienceScripts/DictReplacer.py $INF $IND

