## R script to evaluate sequencing depth in metataxonomics experiments
## custom recursive function to compute the incremental frequency of OTUs
## + code to plot results

library("ggplot2")
library("data.table")

## parameters
basefolder = "/home/filippo/Documents/cremonesi/sheep_f1_methionine/"
inputfile = "Analysis/results/otutable.txt"
outputfile = "Analysis/results/figures/sequencing_depth.png"


#######################################################
## recursive function to do the same loop as above
getIncrementalOtus <- function(tab) {
  
  indx <- list("ind"=which.max(colSums(tab>0)),"max"=max(colSums(tab>0)))
  steps <- indx$max
  
  tab <- tab[tab[,indx$ind]==0,]
  
  if(nrow(tab)>1) {
    
    steps <- c(steps,getIncrementalOtus(tab))
  } else {
    
    res <- ifelse(sum(steps)<nrow(tab),nrow(tab)-sum(steps),0)
    if(res==0) res <- NULL
    steps <- c(steps,res)
  }
  return(steps)
}
#######################################################

## otu table from closed_otupicking (unfiltered)
fname = file.path(basefolder, inputfile)
otu <- fread(fname, header=TRUE, skip=0)
otu <- otu[,-1, with=FALSE]

seqs <- getIncrementalOtus(as.data.frame(otu))
dd <- cbind.data.frame("sample"=seq(1,ncol(otu)),"otus"=cumsum(c(seqs,rep(0,ncol(otu)-length(seqs)))))

p <- ggplot(dd,aes(x=factor(sample),y=otus, group=1)) + geom_line(linewidth=0.75)
p <- p + xlab("N. of samples") + ylab("")
p <- p + scale_x_discrete(breaks=seq(1, 96, 2))
p <- p + theme(text = element_text(size=10),
               axis.text.x = element_text(size=10))
# p

fname = file.path(basefolder, outputfile)
ggsave(filename = fname, plot = p, device = "png")
