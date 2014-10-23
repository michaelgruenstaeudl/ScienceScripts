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
DESCRIPTION="Shell script for R package P2C2M"
AUTHOR="Michael Gruenstaeudl, PhD"
CONTACT="gruenstaeudl.1@osu.edu"
VERSION="2014.10.22.2000"
USAGE="bash <this_script> <abs_path_to_indir> <xml.file> <num.reps> <use.mpi(T/F)>"

################################################################################

# STEP 1: Intro display and exception handling

echo ""
echo -e " ${green}Title: $TITLE | Version: $VERSION | Author: $AUTHOR${nocolor}"
echo -e " ${yellow}Usage: $USAGE${nocolor}"
echo ""

ABS_PATH_TO_INDIR=$1
XML_FILE=$2
NUM_REPS=$3
USE_MPI=$4
# Using Parameter expansion to remove file extension
INFILE=${XML_FILE%.xml*}
#DATE=$(date +%Y-%b-%d)

# Checking number of arguments
if [[ $# != 4 ]]; then
	echo -e " ${red}ERROR: Incorrect number of arguments.${nocolor}"
	exit
fi

for path in $ABS_PATH_TO_INDIR; do
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
if [ ! -f $ABS_PATH_TO_INDIR/$XML_FILE ]; then 
    echo -e " ${red}ERROR: XML file not found.${nocolor}"    
    exit
fi

################################################################################

# STEP 2: Generating R commands for starBeastPPS analysis

echo ""
echo -e " ${blue}~~~ Start of $INFILE ~~~${nocolor}"
echo ""

# Changing input directory
cd $ABS_PATH_TO_INDIR

# Generating R commands
#echo "options(warn=2, error=recover, verbose=TRUE)" > $INFILE.p2c2m.cmd.R
echo "library(P2C2M)" > $INFILE.p2c2m.cmd.R
echo "$INFILE <- p2c2m.complete(path='$ABS_PATH_TO_INDIR', 
                                xml.file='$XML_FILE',
                                descr.stats='GSI,GTP,NDC,RAY',
                                num.reps=$NUM_REPS,
                                use.mpi=$USE_MPI,
                                verbose=TRUE)" >> $INFILE.p2c2m.cmd.R
echo "save($INFILE, file='$INFILE.rda')" >> $INFILE.p2c2m.cmd.R
#echo "warnings()" >> $INFILE.p2c2m.cmd.R
echo "q()" >> $INFILE.p2c2m.cmd.R

# Executing commands in R with command args
Rscript $INFILE.p2c2m.cmd.R

################################################################################

# STEP 3: File hygiene

#echo ""
#echo -e " ${blue}Conducting file hygiene ...${nocolor}"
#echo ""

#tar czf $INFILE.P2C2M.prmt.R.tar.gz $INFILE.P2C2M.prmt.R
#rm $INFILE.P2C2M.prmt.R 
##rm *.P2C2M.phyb

## Check if folder "output" exists on previous level
## Remember: Currently, you are in $ABS_PATH_TO_INDIR
#if [ ! -d ../output ]; then
#  mkdir ../output
#fi
## Separate results from input data
#mv *.P2C2M.* ../output

################################################################################
# EXAMPLE INSTALLATION
#install.packages("P2C2M", dependencies = c("Depends", "Suggests"))

# EXAMPLE EXECUTION
#nohup sh -c '
#for dir in $(ls -d sim* | sort -n -k 1.5); do
#time1=$(date +%Y.%m.%d.%H:%M:%S);
#bash /home/mgruenstaeudl/P2C2MWrapper.sh $(pwd)/$dir/ $dir.xml 100 T;
#time2=$(date +%Y.%m.%d.%H:%M:%S);
#echo "Start time: $time1" >${dir%.xml*}.timelog;
#echo "Stop time:  $time2" >>${dir%.xml*}.timelog;
#done ' &
################################################################################

echo ""
echo -e " ${green}~~~ End of $INFILE ~~~${nocolor}"
echo ""

