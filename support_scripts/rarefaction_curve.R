
# BiocManager::install("dada2")
## Rscript to plot rarefaction curves per sample

library("dada2")
library("vegan")
library("data.table")

## parameters
basefolder = "/home/filippo/Documents/cremonesi/sheep_f1_methionine/"
inputfile = "Analysis/results/otutable.txt"
metadata = "Microbiome-F1-MET-sheep/data/metadata.csv"
outputfile = "Analysis/results/figures/rarefaction_curve.png"

## otu table from closed_otupicking (unfiltered)
fname = file.path(basefolder, inputfile)
otu <- fread(fname, header=TRUE, skip=0)
otu <- otu[,-1, with=FALSE]

## metadata
fname = file.path(basefolder, metadata)
metadata = fread(fname)

## plot and save the rarefaction curves
fname = file.path(basefolder, outputfile)

png(fname, width = 1200, height = 900, res = 250)
rarecurve(t(otu), step=100, col = metadata$MISEQ, lwd=2, xlab = "N. of sequences", ylab="ASVs", label=FALSE)
# and adding a vertical line at the fewest seqs in any sample
abline(v=(min(rowSums(t(otu)))))
dev.off()

