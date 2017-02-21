library(ggtree)
library(svglite) # For improved svg drivers
library(tcltk)


########################
# INFILE / OUTFILE OPS #
########################
# SPECIFYING INFILES
inFile = tk_choose.files(caption='Select .tre file')

# SPECIFYING OUTFILES
# raxml_file <- system.file('~/Desktop/TESTFELD/', 'RAxML_bipartitionsBranchLabels.tmp', package='treeio')
outDir = dirname(inFile)
outFile1 = paste(outDir, '/', 'unrooted_viz.svg', sep='')
outFile2 = paste(outDir, '/', 'rooted_viz.svg', sep='')

# LOAD ALIGNMENTS
tree <- read.raxml(inFile)

# PLOT ML TREES
svglite(outFile1, width=10, height=40, standalone=TRUE)
    ggtree(tree, size=0.75) +
        #geom_label(aes(label=bootstrap), hjust=1, vjust=-0.25) +
        geom_text(aes(label=bootstrap), hjust=1, vjust=-0.25) +
        geom_tiplab() +
        ggtitle('Best ML tree\nModel: GTRGAMMAI, 1000 Bootstrap Replicates\nunrooted')
dev.off()

