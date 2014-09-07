#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2014 Michael Gruenstaeudl"
#email = "gruenstaeudl.1@osu.edu"
#version = "2014.09.05.2000"

## DEBUG Templates
  #cat("\n\n***DEGUB - START***\n")
  #cat("variable")
  #print(variable)
  #cat("\n***DEGUB - END***\n\n")

makeHudsonMSCmd = function(spTree, popSize=5, nAlleles=10) {

  # Description:    Generating gene trees via Hudson's ms   
  # Dependencies:   -
  # InVariables:    spTree:     the species tree; tree must be of class "phylo"
  #                 popSize:    the average population size
  #                 nAlleles:   number of alleles to be generated
  # Notes:          REQUIREMENT: the outgroup MUST be the last tip.label

  library('ape')
  library('adephylo')

  ## Generate merger martix - part 1
    # List the tips of every node
    tipsList = listTips(spTree)
    # Initialize list
    mergers = terminals = c()
    # Loop through all nodes
    for (node in tipsList) {
        # Take the first tip per node
        firstTip = sort(node)[[1]]
        # Take the last tip per node
        lastTip = tail(sort(node), n=1)
        # Ensuring that every i has an existing j to merge with
        # The above problem would occur, for example, if two monophyletic 
        # groups of same rank.
        if (!lastTip %in% terminals) {
            mergers = c(mergers, paste(lastTip, firstTip))
            terminals = c(terminals, lastTip)
        }
        else {
            mergers = c(mergers, paste(firstTip, lastTip))
        }
        if (!firstTip %in% terminals) { terminals = c(terminals, firstTip) }
    }

  ## Generate merger martix - part 2
    # Example: -ej brTms_ij subpop_i subpop_j
    brTms = branching.times(spTree)
    ejList = c()
    for (index in 1:length(brTms)) {
        ejList = c(ejList, paste("-ej", brTms[[index]], mergers[[index]]))
    }

  ## Generate population matrix
    # Generate a list with as many as entries as the tree has tip.labels
    # and fill them with popsize=5
    popList = rep(popSize, length(spTree$tip.label))
    # Give column names to list entries
    names(popList) = paste("pop", 1:length(spTree$tip.label), sep="")
    # Last list entry constitutes outgroup; hence, give it popsize=1
    popList[length(popList)] = 1

  ## Define number and size of populations
    # Example with three populations, each of which has size 5: -I 3 5 5 5
    popInfo = paste("-I", length(spTree$tip.label), paste(popList, collapse=" "))

  ## Define theta
    theta = "-t 1.0"

  ## Assemble the full command
    cmd = paste("ms", sum(popList), nAlleles, theta, popInfo, paste(ejList, collapse=" "), "-T")
    # EXAMPLE OUTPUT: ms <nsam> 10 -t 1.0 -I 4 <pop1> <pop2> <pop3> 1 
    #                 -ej <tau2> 2 1 -ej <tau1> 3 1 -ej <tau3> 4 1 -T

return(cmd)
}
