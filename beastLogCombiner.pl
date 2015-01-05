#!/usr/bin/perl

# developed by Sarah Hird, PhD: shird1@tigers.lsu.edu
# modified/improved by Michael Gruenstaeudl, PhD: gruenstaeudl.1@osu.edu
# version: 14-APR-2014
# This script discards the burn-in of, thins and combines output files from *BEAST
# burnin is the number of steps to be discarded as burn-in
# thinning is the interval at which to subsample MCMC steps for downstream analysis; 
# 10k steps takes quite a while to analyze and may not be necessary. 

use strict;
use Term::ANSIColor;
use warnings;

# Explaining usage
print color("green"), "\nUsage: script.pl <path_to_input> <path_to_output> <burn-in value> <subsampling every>\n\n", color("reset");

my $indirname = $ARGV[0];
my $outdirname = $ARGV[1];
my $burnin = $ARGV[2];
my $thinning = $ARGV[3];

chomp $indirname;
chomp $outdirname;
unless ( -d "$outdirname" )                                                     # make outdir, unless it already exists 
    {mkdir "$outdirname";}

opendir(INDIR, $indirname) or die "can't open dir $indirname: $!";
my @allFiles = grep {!/^\./} grep {!/\.ops/} grep {!/\.xml/} readdir(INDIR);    # grep {!/^\./} means to exclude all files starting with a dot
chomp (@allFiles);
close (INDIR);
my @allFilesSorted = sort(@allFiles);

my $outlogName = "$outdirname/combined.log";
my $treeStatesValue;
my @treeStates= ();
my %current = ();
my $logStatesValue;
my $logCount=0;

foreach my $i (@allFilesSorted)                                                 # for each file in the directory, determine if it is a .trees or .log file, then do something
    {

	if ($i =~ m/.trees/)
        {
		print "  Working on: $i\n";
		@treeStates=();                                                         # clear treeStates array
		my @nameSplit = split(/\./, $i);                                        # split the file name at the periods
		if (exists ($current{"$nameSplit[-2]"}))                                # if the name right before .trees is already in the name hash
            {
			open FILE, "<", "$indirname/$i" or die $!;
			open OUTTREE, ">>", "$outdirname/$current{$nameSplit[-2]}[0]" or die $!;
			my $semicolons = 0;
			while (<FILE>)                                                      #read through the file, get through the header part by counting how many semicolons have been passed
                {
				if ($semicolons < 6)
                    {if ($_ =~ /;/){$semicolons++;}}
				else{push (@treeStates, $_);}                                   # once you've passed the header, push the lines onto the treeStates array
			    }
		    $treeStatesValue = $#treeStates;                                    # count total number of lines 
		    for (my $count = $burnin; $count < $treeStatesValue; $count+=$thinning)
                {print OUTTREE $treeStates[$count];}                             # discard burnin, print correctly thinned lines				
	        }	
		else                                                                    #if file name not in hash, put in hash and do same steps as above to file
            {
			@treeStates=();                                                     #clear treeStates array
			$current{$nameSplit[-2]}[0] = "$nameSplit[-2].combined.trees";
			open FILE, "<", "$indirname/$i" or die $!;
			open OUTTREE, ">>", "$outdirname/$current{$nameSplit[-2]}[0]" or die $!;
			my $semicolons = 0;
			while (<FILE>)
                {
				if ($semicolons < 6)
                    {print OUTTREE "$_";	if ($_ =~ /;/){$semicolons++;}}
				else
                    {push (@treeStates, $_);}
                }
			$treeStatesValue = $#treeStates;
			for (my $count = $burnin; $count < $treeStatesValue; $count+=$thinning)
                {print OUTTREE $treeStates[$count];}	
		    }
	    }
		
	elsif ($i =~ m/.log/)
        {                                                                       #for the log files
		print "  Working on: $i\n";
		open FILE2, "<", "$indirname/$i" or die $!;
		open OUTLOG, ">>", "$outlogName" or die $!;
		my @logStates=();                                                       #clear logStates array
		my $lineCount=0;
		if ($logCount == 0)                                                     #if this is the first log file
            {
			while (<FILE2>)                                                     #print header information to combined log file
                {
                if ($lineCount<3)
                    {print OUTLOG "$_"; $lineCount++;}
				else 
                    {push (@logStates, $_);}
			    }
			$logStatesValue = ($#logStates+1);
			for (my $count2 = $burnin; $count2 < $logStatesValue; $count2+=$thinning)
                {print OUTLOG $logStates[$count2];}
		    $logCount=1;                                                        #switch counter for log file
		    }
		else
            {
			@logStates=();                                                      #clear logStates array
			open FILE2, "<", "$indirname/$i" or die $!;
			my @file2 = <FILE2>;
			$logStatesValue = ($#file2+1);
			for (my $count3 = (3+$burnin); $count3 < $logStatesValue; $count3+=$thinning)
                {print OUTLOG $file2[$count3];}
		    }
	    }

	else
        {print "  Working on: $i -- Ignoring file ...\n";}                      #ignore all other files in directory
	
    }


# Adding an "end;" to the end of the combined trees files
print "\n  Adding an 'end;' to the end of the combined trees files\n";
system("for file in \$(ls $outdirname\*.trees); do echo -e '\nend;\n' >> \$file; done");

# Exiting
print color("green"), "\nDone.\n\n", color("reset");

