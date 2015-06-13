#!/bin/bash

# Global variables
red='\e[1;31m'
blue='\e[1;34m'
green='\e[1;32m'
yellow='\e[0;33m'
purple='\e[0;35m'
cyan='\e[0;36m'
nocolor='\e[0m'

# External functions
# Bash Spinner for Long Running Tasks
# Modified from: http://fitnr.com/showing-a-bash-spinner.html
spinner()
{
    local pid=$1
    local delay=0.75
    local spinstr='|/-\'
    while [ "$(ps a | awk '{print $1}' | grep $pid)" ]; do
        local temp=${spinstr#?}
        printf "${blue} %c  " "$spinstr${nocolor}"
        local spinstr=$temp${spinstr%"$temp"}
        sleep $delay
        printf "${blue}  \b\b\b\b\b\b${nocolor}"
    done
    printf "${blue}    \b\b\b\bDONE${nocolor}"
}

TITLE="cpGenome2partitionedMrBayes.sh"
DESCRIPTION="Shell script to convert annotated cp genome to partitioned MrBayes analysis"
AUTHOR="Michael Gruenstaeudl, PhD"
CONTACT="mi.gruenstaeudl@gmail.com"
VERSION="2015.06.14.2300"
USAGE="bash <this_script> <alignment.nex> <annotations.csv> <path_to_modeltest.jar>"

################################################################################

# STEP 0: Intro display and exception handling

echo ""
echo -e " ${green}Title: $TITLE | Version: $VERSION | Author: $AUTHOR${nocolor}"
echo -e " ${green}Description: $DESCRIPTION${nocolor}"
echo -e " ${yellow}Usage: $USAGE${nocolor}"
echo ""
echo -ne " ${blue} Step 0: Checking files ... ${nocolor}"

ALIGNMENT=$1
ANNOTATIONS=$2
MODELTEST=$3 # /home/michael/binaries/jModelTest_2.1.7/jModelTest.jar

# Using Parameter expansion to remove file extension
INFILE=${ALIGNMENT%.nex*}

# Checking number of arguments
if [[ $# != 3 ]]; then
    echo -e " ${red}ERROR: Incorrect number of arguments.${nocolor}"
    exit
fi

#for path in $ABS_PATH_TO_INDIR; do
## Checking if directories exist
#if [ ! -d $path ]; then 
#    echo -e " ${red}ERROR: Directory "$path" not found.${nocolor}"
#    exit
#fi
## Checking if paths to directories are absolute
#if [[ $path != /* ]]; then
#    echo -e " ${red}ERROR: Path "$path" not absolute.${nocolor}"
#    exit
#fi
# Checking if paths to directories end with forward slash
#if [[ $path != */ ]]; then
#    echo -e " ${red}ERROR: Path "$path" does not end with forward slash.${nocolor}"
#    exit
#fi
#done

# Checking if input files exists
if [ ! -f $ALIGNMENT ]; then 
    echo -e " ${red}ERROR: File not found: $ALIGNMENT${nocolor}"
    exit
fi
if [ ! -f $ANNOTATIONS ]; then 
    echo -e " ${red}ERROR: File not found: $ALIGNMENT${nocolor}"
    exit
fi
if [ ! -f $MODELTEST ]; then 
    echo -e " ${red}ERROR: File not found: $MODELTEST${nocolor}"
    exit
fi

# Generate backup folders
if [ ! -d 01_input ]; then 
    mkdir 01_input
fi

if [ ! -d 02_process ]; then 
    mkdir 02_process
fi

if [ ! -d 03_results ]; then 
    mkdir 03_results
fi

# Make backup of input files
cp $ALIGNMENT 01_input
cp $ANNOTATIONS 01_input

sleep 1 &
spinner $!

################################################################################

# STEP 1: Parse partition information

echo ""
echo -ne " ${blue} Step 1: Parse partition information ... ${nocolor}"

# Remove extraneous lines
cat $ANNOTATIONS | grep -vE "Editing History Deletion|misc_feature|primer_bind|unsure|repeat_region" > $ANNOTATIONS.tmp && mv $ANNOTATIONS.tmp $ANNOTATIONS

# Generating R script
echo "
# Usage: 'Rscript this_script.R input.csv'
#print('Rscript this_script.R input.csv')

# Set R so that passed arguments are read
args = commandArgs(trailingOnly = TRUE)

# Set working directory to the current directory
setwd(getwd())

# GLOBAL VARIABLES
inFn = args[1]

inTable = read.csv(inFn)
# TFLs convert all numbers with commas into numeric format
inTable[,'Minimum'] = suppressWarnings(as.numeric(gsub(',','', inTable[,'Minimum'])))
inTable[,'Maximum'] = suppressWarnings(as.numeric(gsub(',','', inTable[,'Maximum'])))

geneNames = unique(inTable[,'Name'])
outTable = as.data.frame(matrix(nrow=length(geneNames), ncol=3))
rownames(outTable) = geneNames
colnames(outTable) = c('Start', 'End', 'Type')
for (name in geneNames) {
   # The max() in TFL line is correct; don't change it to min()
   outTable[name, 'Start'] = suppressWarnings(max(inTable[which(inTable[,'Name']==name),'Minimum'], na.rm = TRUE))
   outTable[name, 'End'] = suppressWarnings(min(inTable[which(inTable[,'Name']==name),'Maximum'], na.rm = TRUE))
   outTable[name, 'Type'] = toString(unique(inTable[which(inTable[,'Name']==name),'Type'])[1])
}
is.na(outTable[,'Start']) = !is.finite(outTable[,'Start'])
is.na(outTable[,'End']) = !is.finite(outTable[,'End'])
# But not TFL (I don’t know why)
#is.na(outTable[,c('Start','End')]) = !is.finite(outTable[,c('Start','End')])
# Problem: Some of the genes are in forward, others in reverse direction. Hence, the start and stop positions may be incorrect. Hence: [see TFL]
outTable = outTable[which(outTable[,'Start']<outTable[,'End']),]
write.csv(outTable, file=paste(inFn, '.filtered', sep=''))
q()
" > ParsePartitionInfo.R

# Executing R script
Rscript ParsePartitionInfo.R $ANNOTATIONS &
spinner $!

# Add missing column name
sed -i 's/\"\"\,\"Start\"/\"Name\"\,\"Start\"/' $ANNOTATIONS.filtered

# Make backup of this step
cp $ANNOTATIONS.filtered 02_process/

# File hygiene
rm ParsePartitionInfo.R

################################################################################

# STEP 2: Separating partitions of NEXUS file

echo ""
echo -ne " ${blue} Step 2: Separating partitions of NEXUS file ... ${nocolor}"

# Generating a Python script that adds partition info to NEXUS file
echo "
# Usage: 'python2 this_script.py input.nex input.csv.filtered'
#print('python2 this_script.py input.nex input.csv.filtered')

def main():
    # Open nexus alignment file and add partition info
    import csv
    import sys
    with open(sys.argv[1], 'a') as Fhandle_alignm:
        Fhandle_alignm.write('\nbegin sets;\n')
        with open(sys.argv[2]) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                outString = 'charset ' + row['Name'].replace(' ','_').replace('-','_').replace('(','').replace(')','') + ' = ' + row['Start'] + '-' + row['End'] + ';\n'
                Fhandle_alignm.write(outString)
        Fhandle_alignm.write('end;\n')

    # Open nexus alignment file (which includes partition info) and split by partition
    from Bio.Nexus import Nexus
    aln = Nexus.Nexus()
    aln.read(sys.argv[1])
    aln.write_nexus_data_partitions(filename='partition', charpartition=aln.charsets)

main()
exit()
" > SplitPartitions.py

# Executing Python script (and immediately thereafter the Bash spinner)
python2 SplitPartitions.py $ALIGNMENT $ANNOTATIONS.filtered &
spinner $!

# Ordering the partitions by adding consecutive numbers
n=1;
for i in partition*; do 
rename partition partition$(printf %04d $n) $i; n=$((n+1));
done

# Make backup of this step
cp partition* 02_process/

# File hygiene
rm SplitPartitions.py

################################################################################

# STEP 3: Conducting modeltesting via jModelTest2

echo ""
echo -ne " ${blue} Step 3: Conducting modeltesting via jModelTest2 ... ${nocolor}"

# Generating a Python script that adds partition info to NEXUS file
for ds in $(ls partition*); do 
java -jar $MODELTEST -d $ds -g 4 -i -f -AIC -o $ds.bestModel 1>$ds.bestModel.log 2>$ds.bestModel.err &; 
spinner $!;
done

# Extract model information from files
for ds in $(ls *.bestModel); do echo $ds >> model_overview.txt; cat $ds | grep -A1 ' Model selected:' | tail -n1 >> model_overview.txt; done

# Make backup of this step
cp *.bestModel 02_process/
mv *.bestModel.log 02_process/
mv *.bestModel.err 02_process/
cp model_overview.txt 02_process/

################################################################################

# STEP 4: Setting up a partitioned MrBayes analysis

echo ""
echo -ne " ${blue} Step 4: Setting up a partitioned MrBayes analysis ... ${nocolor}"

# Generating a Python script to set up a partitioned nexus file
echo "
# Usage: 'python2 this_script.py'
#print('python2 this_script.py')

def main():
    from Bio.Nexus import Nexus
    import glob
    file_list = glob.glob('partition*')
    file_list.sort() # Important step
    nexi = [(fname, Nexus.Nexus(fname)) for fname in file_list] 
    combined = Nexus.combine(nexi)
    combined.write_nexus_data(filename=open('combined.nex', 'w'))

main()
exit()
" > CombinePartitions.py

# Executing Python script (and immediately thereafter the Bash spinner)
python2 CombinePartitions.py &
spinner $!

# Setting up lset specifications
sed -i 's/\.nex \=//g' combined.nex
sed -i ':a;N;$!ba;s/\.nex\.bestModel\n/\)/g' model_overview.txt
sed -i ':a;N;$!ba;s/\.nex\.bestModel/\)/g' model_overview.txt

sed -i 's/Model \= F81\+I\+G/nst\=1 rates\=invgamma\;/' model_overview.txt
sed -i 's/Model \= JC\+I\+G/nst\=1 rates\=invgamma\;/' model_overview.txt
sed -i 's/Model \= K80\+I\+G/nst\=2 rates\=invgamma\;/' model_overview.txt
sed -i 's/Model \= HKY\+I\+G/nst\=2 rates\=invgamma\;/' model_overview.txt
sed -i 's/Model \= SYM\+I\+G/nst\=6 rates\=invgamma\;/' model_overview.txt
sed -i 's/Model \= GTR\+I\+G/nst\=6 rates\=invgamma\;/' model_overview.txt

sed -i 's/Model \= F81\+G/nst\=1 rates\=gamma\;/' model_overview.txt
sed -i 's/Model \= JC\+G/nst\=1 rates\=gamma\;/' model_overview.txt
sed -i 's/Model \= K80\+G/nst\=2 rates\=gamma\;/' model_overview.txt
sed -i 's/Model \= HKY\+G/nst\=2 rates\=gamma\;/' model_overview.txt
sed -i 's/Model \= SYM\+G/nst\=6 rates\=gamma\;/' model_overview.txt
sed -i 's/Model \= GTR\+G/nst\=6 rates\=gamma\;/' model_overview.txt

sed -i 's/Model \= F81\+I/nst\=1 rates\=propinv\;/' model_overview.txt
sed -i 's/Model \= JC\+I/nst\=1 rates\=propinv\;/' model_overview.txt
sed -i 's/Model \= K80\+I/nst\=2 rates\=propinv\;/' model_overview.txt
sed -i 's/Model \= HKY\+I/nst\=2 rates\=propinv\;/' model_overview.txt
sed -i 's/Model \= SYM\+I/nst\=6 rates\=propinv\;/' model_overview.txt
sed -i 's/Model \= GTR\+I/nst\=6 rates\=propinv\;/' model_overview.txt

sed -i 's/Model \= F81/nst\=1\;/' model_overview.txt
sed -i 's/Model \= JC/nst\=1\;/' model_overview.txt
sed -i 's/Model \= K80/nst\=2\;/' model_overview.txt
sed -i 's/Model \= HKY/nst\=2\;/' model_overview.txt
sed -i 's/Model \= SYM/nst\=6\;/' model_overview.txt
sed -i 's/Model \= GTR/nst\=6\;/' model_overview.txt

for keyw in $(cat combined.nex | grep charset | awk '{print $2}'); do 
linenum=$(cat combined.nex | grep charset | awk '/'$keyw'/{print NR; exit}');
sed -i "s/$keyw/$linenum/" model_overview.txt; # Double quotes are critical
done

awk '{print "lset apply=(" $0}' model_overview.txt > model_overview.txt.tmp && mv model_overview.txt.tmp model_overview.txt

# Setting up MrBayes nexus
CHARSET_STRING="$(cat combined.nex | grep charset | awk '{print $1" "$2" = "$3}')"
cat combined.nex | grep -v charset | grep -v charpartition | sed '$d' | sed '$d' >> combined.nex.tmp && mv combined.nex.tmp combined.nex
echo 'BEGIN MRBAYES;' >> combined.nex
echo $CHARSET_STRING >> combined.nex
sed -i ':a;N;$!ba;s/\; charset/\;\ncharset/g' combined.nex
echo -n 'partition combined = ' >> combined.nex
echo -n $(cat combined.nex | grep charset | wc -l) >> combined.nex
echo -n ': ' >> combined.nex
cat combined.nex | grep charset | awk '{print $2}' | awk '{ 
ORS=", "; print; }' >> combined.nex
echo ';' >> combined.nex
sed -i 's/\, \;/\;/g' combined.nex
echo 'set partition = combined;' >> combined.nex

cat model_overview.txt >> combined.nex
echo 'mcmcp ngen=1000000 temp=0.1 samplefreq=10000; mcmc;' >> combined.nex
echo 'END; QUIT;' >> combined.nex

# Make backup of this step
cp combined.nex 02_process/

# File hygiene
rm CombinePartitions.py

################################################################################

echo ""
echo -e " ${green}~~~ End of $INFILE ~~~${nocolor}"
echo ""
