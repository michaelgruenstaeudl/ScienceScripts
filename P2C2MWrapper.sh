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
VERSION="2014.09.20.1400"
USAGE="bash <this_script> <abs_path_to_indir> <abs_path_to_scrptdir> <abs_path_to_MS> <xml_infile> <nreps_flag(integer)> <debug_flag(T/F)>"

################################################################################

# STEP 1: Intro display and exception handling

echo ""
echo -e " ${green}Title: $TITLE | Version: $VERSION | Author: $AUTHOR${nocolor}"
echo -e " ${yellow}Usage: $USAGE${nocolor}"
echo ""

ABS_PATH_TO_INDIR=$1
ABS_PATH_TO_SCRIPTDIR=$2
ABS_PATH_TO_MS=$3
XML_INFILENAME=$4
NREPS_FLAG=$5
DEBUG_FLAG=$6

# Using Parameter expansion to remove file extension
INFILE=${XML_INFILENAME%.xml*}
#DATE=$(date +%Y-%b-%d)

# Checking number of arguments
if [[ $# != 6 ]]; then
	echo -e " ${red}ERROR: Incorrect number of arguments.${nocolor}"
	exit
fi

for path in $ABS_PATH_TO_INDIR $ABS_PATH_TO_SCRIPTDIR $ABS_PATH_TO_MS; do
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
echo -e " ${blue}~~~ Analyzing $INFILE ~~~${nocolor}"
echo ""

# Changing input directory
cd $ABS_PATH_TO_INDIR

# Generating R commands
echo "source('$ABS_PATH_TO_SCRIPTDIR/wrapper.R')" > $INFILE.P2C2M.Rcmds
echo "wrapper.init( '$ABS_PATH_TO_INDIR', 
                    '$XML_INFILENAME',
                    '$ABS_PATH_TO_SCRIPTDIR',
                    '$ABS_PATH_TO_MS',
                    '$NREPS_FLAG',
                    '$DEBUG_FLAG'   )" >> $INFILE.P2C2M.Rcmds
echo "warnings()" >> $INFILE.P2C2M.Rcmds
echo "q()" >> $INFILE.P2C2M.Rcmds

# Executing commands in R with command args
Rscript $INFILE.P2C2M.Rcmds

# Correction of .P2C2M.rslts.csv
sed -i 's/"GTP\[2\]"/,"GTP\[2\]"/g' $INFILE.P2C2M.rslts.csv

################################################################################

# STEP 3: File hygiene

echo ""
echo -e " ${blue}Conducting file hygiene ...${nocolor}"
echo ""

tar czf $INFILE.P2C2M.prmt.R.tar.gz $INFILE.P2C2M.prmt.R
rm $INFILE.P2C2M.prmt.R *.P2C2M.phyb

# Check if folder "output" exists on previous level
# Remember: Currently, you are in $ABS_PATH_TO_INDIR
if [ ! -d ../output ]; then
  mkdir ../output
fi
# Separate results from input data
mv *.P2C2M.* ../output

################################################################################

echo ""
echo -e " ${green}~~~ End of Script ~~~${nocolor}"
echo ""

