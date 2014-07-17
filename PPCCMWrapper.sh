#!/bin/bash

# Global variables
red='\e[1;31m'
blue='\e[1;34m'
green='\e[1;32m'
yellow='\e[0;33m'
purple='\e[0;35m'
cyan='\e[0;36m'
nocolor='\e[0m'

TITLE="PPCCMWrapper.sh"
DESCRIPTION="Shell script that wraps the commands necessary to perform a posterior predictive checking of the coalescent model"
AUTHOR="Michael Gruenstaeudl, PhD"
CONTACT="gruenstaeudl.1@osu.edu"
VERSION="2014.07.17.1700"
USAGE="bash <this_script> <abs_path_to_indir> <abs_path_to_outdir> <abs_path_to_MS> <xml_infile> <debug_flag(0/1)>"

################################################################################

# STEP 1: Intro display and exception handling

echo ""
echo -e " ${green}Title: $TITLE | Version: $VERSION | Author: $AUTHOR${nocolor}"
echo -e " ${yellow}Usage: $USAGE${nocolor}"
echo ""

ABS_PATH_TO_INDIR=$1
ABS_PATH_TO_OUTDIR=$2
ABS_PATH_TO_MS=$3
XML_INFILENAME=$4
DEBUG_FLAG=$5

# Using Parameter expansion to remove file extension
INFILE=${XML_INFILENAME%.xml*}
#DATE=$(date +%Y-%b-%d)

# Checking number of arguments
if [[ $# != 5 ]]; then
	echo -e " ${red}ERROR: Incorrect number of arguments.${nocolor}"
	exit
fi

for path in $ABS_PATH_TO_INDIR $ABS_PATH_TO_OUTDIR $ABS_PATH_TO_MS; do
# Checking if directories exist
if [ ! -d $path ]; then 
    echo -e " ${red}ERROR: Directory "$path" not found.${nocolor}"    
    exit
fi
# Checking if paths to directories are absolute
if [[ $path != /* ]]; then
    echo -e " ${red}ERROR: Path "$path" not absolute.${nocolor}"
    exit
fi
# Checking if paths to directories end with forward slash
if [[ $path != */ ]]; then
    echo -e " ${red}ERROR: Path "$path" does not end with forward slash.${nocolor}"
    exit
fi
done

# Checking if input file exists
if [ ! -f $ABS_PATH_TO_INDIR/$XML_INFILENAME ]; then 
    echo -e " ${red}ERROR: XML file not found.${nocolor}"    
    exit
fi

################################################################################

# STEP 2: Generating R commands for starBeastPPS analysis

echo " Analyzing $INFILE ..."

# Changing input directory
cd $ABS_PATH_TO_INDIR

# Generating R commands
echo "source('/home/michael/git/osu_eeob_gruenstaeudl/starbeast_Michael/wrapper.R')" > Rcmds.$INFILE.R
echo "wrapper.go('$ABS_PATH_TO_INDIR','$XML_INFILENAME','$ABS_PATH_TO_MS')" >> Rcmds.$INFILE.R
#echo "warnings()" >> Rcmds.$INFILE.R
echo "q()" >> Rcmds.$INFILE.R

# Executing commands in R with command args
Rscript Rcmds.$INFILE.R $DEBUG_FLAG

################################################################################

# STEP 3: File hygiene
mv $INFILE.PRMT $INFILE.RSLT $ABS_PATH_TO_OUTDIR/
rm Rcmds.$INFILE.R
cd ~
#tar czf ${ABS_PATH_TO_INDIR:0:-1}.tar.gz $ABS_PATH_TO_INDIR

################################################################################

echo ""
echo -e " ${green}--- End of Script ---${nocolor}"
echo ""

