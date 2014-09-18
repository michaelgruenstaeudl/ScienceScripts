#!/bin/bash

# Global variables
red='\e[1;31m'
blue='\e[1;34m'
green='\e[1;32m'
yellow='\e[0;33m'
purple='\e[0;35m'
cyan='\e[0;36m'
nocolor='\e[0m'

TITLE="P2C2MWrapper.sh"
DESCRIPTION="Shell script that wraps the commands necessary to perform a posterior predictive checking of the coalescent model"
AUTHOR="Michael Gruenstaeudl, PhD"
CONTACT="gruenstaeudl.1@osu.edu"
VERSION="2014.09.17.2300"
USAGE="bash <this_script> <abs_path_to_indir> <abs_path_to_outdir> <abs_path_to_scrptdir> <abs_path_to_MS> <xml_infile> <nreps_flag(integer)> <debug_flag(T/F)>"

################################################################################

# STEP 1: Intro display and exception handling

echo ""
echo -e " ${green}Title: $TITLE | Version: $VERSION | Author: $AUTHOR${nocolor}"
echo -e " ${yellow}Usage: $USAGE${nocolor}"
echo ""

ABS_PATH_TO_INDIR=$1
ABS_PATH_TO_OUTDIR=$2
ABS_PATH_TO_SCRIPTDIR=$3
ABS_PATH_TO_MS=$4
XML_INFILENAME=$5
NREPS_FLAG=$6
DEBUG_FLAG=$7

# Using Parameter expansion to remove file extension
INFILE=${XML_INFILENAME%.xml*}
#DATE=$(date +%Y-%b-%d)

# Checking number of arguments
if [[ $# != 7 ]]; then
	echo -e " ${red}ERROR: Incorrect number of arguments.${nocolor}"
	exit
fi

for path in $ABS_PATH_TO_INDIR $ABS_PATH_TO_OUTDIR $ABS_PATH_TO_SCRIPTDIR $ABS_PATH_TO_MS; do
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

echo ""
echo -e " ${blue}Analyzing $INFILE ...${nocolor}"
echo ""

# Changing input directory
cd $ABS_PATH_TO_INDIR

# Generating R commands
echo "source('$ABS_PATH_TO_SCRIPTDIR/wrapper.R')" > Rcmds.$INFILE.R
echo "wrapper.go('$ABS_PATH_TO_INDIR','$XML_INFILENAME', '$ABS_PATH_TO_SCRIPTDIR', '$ABS_PATH_TO_MS', '$NREPS_FLAG', '$DEBUG_FLAG')" >> Rcmds.$INFILE.R
echo "warnings()" >> Rcmds.$INFILE.R
echo "q()" >> Rcmds.$INFILE.R

# Executing commands in R with command args
Rscript Rcmds.$INFILE.R

################################################################################

# STEP 3: File hygiene

echo ""
echo -e " ${blue}Conducting file hygiene ...${nocolor}"
echo ""

tar czf $INFILE.PRMT.tar.gz $INFILE.PRMT
rm $INFILE.PRMT
mv $INFILE.PRMT.tar.gz $INFILE.RSLT $INFILE.*.R $ABS_PATH_TO_OUTDIR/

cd ~

################################################################################

echo ""
echo -e " ${green}--- End of Script ---${nocolor}"
echo ""

