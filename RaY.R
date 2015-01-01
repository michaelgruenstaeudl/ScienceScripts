#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2014 Michael Gruenstaeudl"
#email = "gruenstaeudl.1@osu.edu"
#version = "2014.06.16.1700"

###################
# DEBUG Templates #
###################
#cat("\n\n***DEGUB - START***\n")
#cat("variable")
#print(variable)
#cat("\n***DEGUB - END***\n\n")

#######################
# Function "RaY.main" #
#######################
# Description:      Administrating/organizing the Ranalla&Yang calculations
# Subsidiaries:     RaY.calc
# Input:            gtFn = gene tree file name
#                   stFn = species tree file name
#                   pySc = absolute path to Python script "parse_sptree_beast_to_phybase.py"
# Notes:            -

RaY.main = function (gtFn, stFn, pySc="/home/michael/git/CustomScripts/beastTree2phybaseTree.py")
{
#DEBUGLINE cat("\n    > RaY.main\n")

    library('ape')
    library('apTreeshape')
    library('phybase')
    library('phyloch')
    library('xtermStyle')
    source('/home/michael/git/CustomScripts/GeneralStringOperations.R')

    # TFL is not needed, as current wd is set automatically
    #setwd(getwd())

#1. Loading of gene trees
#    gtrees = tryCatch(
#        {tree=ape::read.nexus(gtFn); return(tree)}, 
#        error=function(e) {tree=phyloch::read.beast(gtFn); return(tree)}
#                    )
    # TFL is just a temporary fix; atcually there should be a tryCatch()
    gtrees = ape::read.nexus(gtFn)

#2.a Loading of species trees
    # TFL executes the Python script on the species tree
    system(paste('python2',pySc,'-i',stFn))
    stFN.new = rmext(stFn)
    streeData = phybase::read.tree.string(paste(stFN.new,'.phyb',sep=''))
    strees = streeData$tree
    streeNames = streeData$names

#2.b Special parsing of species tree names to remove "+E" characters
    newStreeNames = c()
    for (i in 1:length(streeNames))
        {
        if ("+" %in% strsplit(streeNames[i],"")[[1]])
            {
            cat("\n",style("WARNING: Taxon names of species tree had to be adjusted due to presence of E-values in branch lengths.",fg="orange"),"\n",sep="")
            # TFL provides index number of keyword in name; [1] is important, because 'which' returns multiple cases
            pos = which(match(strsplit(streeNames[i],"")[[1]],"+")==1)[1]
            streeNames[i] = paste(unlist(strsplit(streeNames[i],""))[1:pos-1],collapse="")
            }
        if ("-" %in% strsplit(streeNames[i],"")[[1]])
            {
            cat("\n",style("WARNING: Taxon names of species tree had to be adjusted due to presence of E-values in branch lengths.",fg="orange"),"\n",sep="")
            pos = which(match(strsplit(streeNames[i],"")[[1]],"-")==1)[1]
            streeNames[i] = paste(unlist(strsplit(streeNames[i],""))[1:pos-1],collapse="")
            }
        }

#3. Decision tree if multiple gene trees and/or multiple species trees

    # Due to the tree reading method, the species trees will always be read in 
    # as a list, even if there is only a single species tree in the file.
    # The gene trees, by contrast, may be read in as a string (if only a single 
    # tree exists in file) or as a list, because ape::read.nexus specifies it as such.

    if (length(gtrees)==1)
        {
        if (length(strees)==1)
            {
            cat("\n",style("DETECTED: Single gene tree and single species tree.",fg="green"),"\n",sep="")
            streeInfo = list()
            streeInfo[[1]] = streeNames
            streeInfo[[2]] = strees
            value = RaY.calc(gtrees, streeInfo)
            results = c(i,"GTsingle-STsingle",value)
            }
        outD = as.data.frame(results)
        colnames(outD)=c("GTnum","compar","LLval")
        if (length(strees)>1)
            {
            cat("\n",style("ERROR: Only a single gene tree, but multiple species trees detected. Are you sure?",fg="red"),"\n",sep="")
            break()
            }
        }

    if (length(gtrees)>1)
        {
        if (length(strees)==1)
            {
            cat("\n",style("DETECTED: Multiple gene trees, but only a single species tree.",fg="green"),"\n",sep="")
            results = c()
            for (i in 1:length(gtrees))
                {
                streeInfo = list()
                streeInfo[[1]] = streeNames
                streeInfo[[2]] = strees
                value = RaY.calc(gtrees[[i]], streeInfo)
                comp = paste("GT",i,"-","STsingle",sep="")
                results = rbind(results, c(i,comp,value))
                }
            outD = as.data.frame(results)
            colnames(outD)=c("GTnum","compar","LLval")
            }
        if (length(strees)>1)
            {
            if (length(gtrees)!=length(strees))
                {
                cat("\n",style("ERROR: Multiple gene trees and multiple species trees, but unequal tree numbers.",fg="red"),"\n",sep="")
                break()
                }
            if (length(gtrees)==length(strees))
                {
                cat("\n",style("DETECTED: Multiple gene trees, multiple species trees, equal tree numbers.",fg="green"),"\n",sep="")
                results = c()
                for (i in 1:length(gtrees))
                    {
                    streeInfo = list()
                    streeInfo[[1]] = streeNames
                    streeInfo[[2]] = strees[i]
                    value = RaY.calc(gtrees[[i]], streeInfo)
                    comp = paste("GT",i,"-","ST",i,sep="")
                    results = rbind(results, c(i,comp,value))
                    }
                outD = as.data.frame(results)
                colnames(outD)=c("GTnum","compar","LLval")
                }
            }
        }
    
    return(outD)
}


#######################
# Function "RaY.calc" #
#######################
# Description:      Conducting the Ranalla&Yang calculation
# Subsidiaries:     -
# Input:            gtFn = gene tree file name
#                   stFn = species tree file name
# Notes:            -

RaY.calc = function (gtree, streeInfo)
{
#1. Loading of gene tree attributes
    gtreeTreeshape = apTreeshape::as.treeshape(gtree)
    gtreeTaxa = gtreeTreeshape$names
    gtreeString = ape::write.tree(gtree)

#2. Loading of species tree attributes
    streeTaxa = streeInfo[[1]]
    streeString = streeInfo[[2]]

#3. Decision structure regarding id labels
#    Explanation: There are two main species-taxon binding options
#    *Opt.1*     *Opt.2*
#    TAX SP      TAX   SP
#    a1  pop1    aaa1  aaa
#    a2  pop1    aaa2  aaa
#    b3  pop2    bbb1  bbb
#    c4  pop3    ccc1  ccc

    idsTaxaPre = unlist(lapply(gtreeTaxa,substring,first=1,last=1))
    idsTaxa = unique(idsTaxaPre)
    idsSpecies = unique(unlist(lapply(streeTaxa,substring,first=1,last=1)))
    spLen = length(unlist(strsplit(streeTaxa[1],"")))
    # Switching to option 2 if number of unique first SP letters <= 2
    if (length(idsSpecies) <= 2)
        {idsSpecies = unique(unlist(lapply(streeTaxa,substring,first=spLen,last=spLen)))}

    # Confirming if number of unique (!) gene tree OTUs equal to number of unique species
    if (length(idsTaxa) != length(idsSpecies))
        {print("ERROR in Function Ranalla&Yang: N gene tree OTUs unequal N species!")
        break}

#4. Generating gene-species dictionary

    counter = list()
    dict = vector(mode="list")
    for (i in 1:length(idsTaxa))
        {
        # FL means "if element not yet in counter list"
        if (! is.null(counter[i]))
            {
            counter = c(counter,i)
            # The resulting dictionary has the gene tree OTUs as keys 
            # and the species as values
            dict[[idsTaxa[[i]]]] = idsSpecies[[i]]
            }
        }

#5. Generating taxon-species binding matrix
    # In TFLs it is important to use idsTaxaPre instead of idsTaxa
    bindings = matrix(0,nrow=length(idsSpecies),ncol=length(idsTaxaPre))
    rownames(bindings) = idsSpecies
    colnames(bindings) = idsTaxaPre
    for (i in 1:length(dict))
        {
        # TFL returns all those columns by index number that are named 'names(dict[i]))'
        # TFL is important, because the code allows multiple columns to have the same label
        collbls = which(colnames(bindings)==names(dict[i]))
        bindings[dict[[i]],collbls] = 1
        }
    # TFLs reassign the correct row and column names to bindings
    rownames(bindings) = streeTaxa
    colnames(bindings) = gtreeTaxa

#6. Performing results calc.
    result = loglikeSP(gtreeString,streeString,gtreeTaxa,streeTaxa,bindings)
    return(result)
}
