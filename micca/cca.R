
library("plyr")
library("broom")
library("dplyr")
library("vegan")
library("caret")
library("ggrepel")
library("tidyverse")
library("tidymodels")
library("data.table")
library("gghighlight")

# library("RANN")
# library("scales")
# library("ggplot2")
# library("madbito")
# library("reshape2")


## PARAMETERS
args = commandArgs(trailingOnly=TRUE)
if (length(args) >= 1) {
  
  #loading the parameters
  if (file_ext(args[1]) %in% c("r","R")) {
    
    source(args[1])
    # source("Analysis/hrr/config.R")
  } else {
    
    load(args[1])
  }
  
} else {
  #this is the default configuration, used for development and debug
  writeLines('Using default config')
  
  #this dataframe should be always present in config files, and declared
  #as follows
  config = NULL
  config = rbind(config, data.frame(
    #base_folder = '~/Documents/SMARTER/Analysis/hrr/',
    #genotypes = "Analysis/hrr/goat_thin.ped",
    repo = "Documents/cremonesi/metabarcoding",
    prjfolder = "Documents/cremonesi/vitelli_giulia_sala_2025",
    analysis_folder = "Analysis",
    conf_file = "Config/mapping_file.csv",
    omics_file = "Analysis/otu_norm_CSS.csv",
    suffix = "calves_colostrum",
    cols_to_keep = "2:113", ## position of omics cols to keep (numeric columns): start end position (separated by :,.-_)
    normalize = TRUE, ## should we normalise omics data?
    nfactors = 1, ## n. of design variables (e.g. treatment and timpoint --> nfactors = 2)
    min_tot_n = 10,
    min_sample = 2,
    project = "", ##! use only for subsetting
    cov_factor = "timepoint",
    cov_continuous = "weight",
    sample_column = "sample_id",
    sample_prefix = "sample-",
    grouping_variable = "treatment",
    base_treatment = 5, ## reference level within timepoint (e.g. control)
    base_timepoint = 1, ## reference level within treatment (e.g. T0)
    force_overwrite = FALSE
  ))
}


HOME <- Sys.getenv("HOME")
repo = file.path(HOME, config$repo)
prjfolder = file.path(HOME, config$prjfolder)
outdir = file.path(prjfolder,config$analysis_folder)

fname = paste("cca.config_",config$suffix,".RData", sep="")
fname = file.path(outdir, fname)
save(config, file = fname)

## METADATA
fname = file.path(HOME, config$prjfolder, config$conf_file)
metadata <- fread(fname)
if(config$sample_prefix != "") {
  
  temp = paste(config$sample_prefix, metadata$sample_id, sep="")
  metadata <- mutate(metadata, !!config$sample_column := temp)
}

### OMICS DATA
fname = file.path(HOME, config$prjfolder, config$omics_file)
omics <- fread(fname)
feature_ids = select(omics ,1) |> pull() ## ordered vector of omics features

## select columns
temp = unlist(strsplit(config$cols_to_keep, ":|,|\\.|-|_"))
position = as.numeric(temp[1]):as.numeric(temp[2])
omics <- omics %>% select(position)

## filter and normalize
omics <- omics[rowSums(omics) > config$min_tot_n,] #keep only small RNA with more than 10 counts

if (config$normalize) {
  
  std = sapply(omics, sd)
  omics <- (omics-colMeans(omics))
  summary(colMeans(omics))
  omics <- sweep(omics, 2, std, FUN="/")
  summary(sapply(omics, mean)) 
}

ids = names(omics)
M <- t(omics)

## DISTANCE MATRIX
D <- dist(M, method="euclidean", diag = TRUE)
D <- data.matrix(D)

heatmap(D,col=heat.colors(75))

## MULTIDIMENSIONAL SCALING
## eigenvalues
A <- cmdscale(D, eig=TRUE, k=3)

###########
## mds plot
###########
mdsD <- as.data.frame(A$points)
names(mdsD) <- c("dim1","dim2","dim3")

mdsD$id = rownames(mdsD)
mdsD <- relocate(mdsD, "id")

temp <- metadata |> 
  select(!!config$sample_column, 
         !!config$grouping_variable, 
         !!config$cov_continuous, 
         !!config$cov_factor)

left_col <- "id"
right_col = config$sample_column
join_vars = setNames(right_col, left_col)

mdsD <- mdsD |> 
  inner_join(temp, by = join_vars)

mdsD <- mdsD |>
  mutate(!!sym(config$cov_factor) := as.factor(.data[[config$cov_factor]]),
         !!sym(config$grouping_variable) := as.factor(.data[[config$grouping_variable]]))


### HULLS
find_hull <- function(mdsD) mdsD[chull(mdsD$dim1, mdsD$dim2), ]
hulls <- ddply(mdsD, "timepoint", find_hull)

hulls <- mdsD %>%
  group_by(timepoint) %>%
  group_modify(~ find_hull(.x)) %>%
  ungroup()

ph <- ggplot(mdsD,aes(x=dim1,y=dim2))
ph <- ph + geom_point(aes(colour=treatment),size=4)
ph <- ph + geom_polygon(data=hulls, alpha=.2, aes(fill=timepoint))
# ph <- ph + scale_fill_discrete("Status",labels=c("euthyroid","hyperthyroid"))
ph <- ph + labs(colour = "treatment") #+ ggtitle("MDS of miRNA+proteomics-based distances")
ph <- ph + theme(
  legend.title=element_text(size=22),
  legend.text=element_text(size=18),
  axis.title=element_text(size=17),
  axis.text=element_text(size=14))
ph

### CCA
M <- mm[,c(2:4)]
H <- mm[,c(6:10)]
preProcValues <- preProcess(H, method = c("knnImpute"))
H <- predict(preProcValues, H)
# M <- t(apply(M,1,rescale))

simX <- as.matrix(cMP[,-1])

vec <- sample(ncol(simX), 20)
X <- simX[,vec]

cca1 <- cca(X ~ TRAB+TSH+FT4, data=H, scale=TRUE)

V <- as.data.frame(cca1$CCA$wa)
V$id <- row.names(V)
V <- merge(V,mm[,c(1,5,6,10)],by="id")

B <- as.data.frame(cca1$CCA$biplot)

pCCA <- ggplot(V,aes(CCA1,CCA2)) + geom_point(aes(colour=status,size=3))
pCCA <- pCCA + geom_segment(data = B, aes(xend = B[ ,"CCA1"], yend=B[ ,"CCA2"]),
                            x=0, y=0, colour="black",
                            arrow=arrow(angle=25, length=unit(0.5, "cm")))
pCCA <- pCCA + geom_text(data=B, aes(x=B[ ,"CCA1"], y=B[ ,"CCA2"], label=row.names(B)), 
                         size=6, vjust=1.5, colour="black")
pCCA <- pCCA + theme(strip.text = element_text(size=20),
                     axis.title=element_text(size=17),
                     axis.text=element_text(size=14))
pCCA <- pCCA + theme(legend.position="none")


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

 
multiplot(pCCA,ph,cols=2)







