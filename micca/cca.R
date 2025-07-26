
library("plyr")
library("broom")
library("dplyr")
library("vegan")
library("caret")
library("ggpubr")
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
    normalize = FALSE, ## should we normalise omics data?
    distance = "minkowski",
    min_tot_n = 50,
    min_sample = 2,
    project = "", ##! use only for subsetting
    cov_factor = "fecal_score",
    cov_continuous = "timepoint,weight",
    sample_column = "sample_id",
    sample_prefix = "sample-",
    grouping_variable = "treatment",
    base_treatment = 5, ## reference level within timepoint (e.g. control)
    base_timepoint = 1, ## reference level within treatment (e.g. T0)
    force_overwrite = FALSE
  ))
}

writeLines(" - set up")
HOME <- Sys.getenv("HOME")
repo = file.path(HOME, config$repo)
prjfolder = file.path(HOME, config$prjfolder)
outdir = file.path(prjfolder,config$analysis_folder)

writeLines(" - saving configs")
fname = paste("cca.config_",config$suffix,".RData", sep="")
fname = file.path(outdir, fname)
save(config, file = fname)

## METADATA
writeLines(" - read metadata")
fname = file.path(HOME, config$prjfolder, config$conf_file)
metadata <- fread(fname)
if(config$sample_prefix != "") {
  
  temp = paste(config$sample_prefix, metadata$sample_id, sep="")
  metadata <- mutate(metadata, !!config$sample_column := temp)
}

### OMICS DATA
writeLines(" - read omics data")
fname = file.path(HOME, config$prjfolder, config$omics_file)
omics <- fread(fname)
feature_ids = select(omics ,1) |> pull() ## ordered vector of omics features

writeLines(" - data preprocessing")
## select columns
temp = unlist(strsplit(config$cols_to_keep, ":|,|\\.|-|_"))
position = as.numeric(temp[1]):as.numeric(temp[2])
omics <- omics %>% select(position)

## filter
omics <- omics[rowSums(omics) > config$min_tot_n,] #keep only small RNA with more than 10 counts

## normalize
new_scale <- function(x) {x/sum(x)}
if (config$normalize) {
  
  # std = sapply(omics, sd)
  # omics <- (omics-colMeans(omics))
  # summary(colMeans(omics))
  # omics <- sweep(omics, 2, std, FUN="/")
  # summary(sapply(omics, mean)) 
  
  omics <- omics |>
    mutate(across(where(is.numeric), scale))
  
  # omics <- omics |>
  #   mutate(across(where(is.numeric), new_scale))
}

ids = names(omics)
M <- t(omics)

## DISTANCE MATRIX
writeLines(" - calculating distances")
D <- dist(M, method=config$distance, diag = TRUE)
D <- data.matrix(D)

fname = file.path(outdir, "heatmap_distance.png")
png(fname)
heatmap(D, col=heat.colors(75))
dev.off()

## MULTIDIMENSIONAL SCALING
writeLines(" - Multidimensional scaling of the distance matrix")
## eigenvalues
A <- cmdscale(D, eig=TRUE, k=3)

###########
## mds plot
###########
mdsD <- as.data.frame(A$points)
names(mdsD) <- c("dim1","dim2","dim3")

mdsD$id = rownames(mdsD)
mdsD <- relocate(mdsD, "id")

if(grepl(pattern = ",", x = config$cov_continuous)) {
  
  continuous_covs = unlist(strsplit(config$cov_continuous,split = ","))
} else continuous_covs = config$cov_continuous

## impute categorical variable
fct_x = select(metadata, !!config$cov_factor, all_of(continuous_covs))
preProcValues <- preProcess(fct_x, method = c("medianImpute"))
imputed_fecal_score = predict(preProcValues, select(fct_x, !!config$cov_factor)) |> pull()
fct_x <- mutate(fct_x, !!config$cov_factor := imputed_fecal_score)

temp <- metadata |> 
  select(!!config$sample_column, 
         !!config$grouping_variable)

temp <- bind_cols(temp, fct_x)

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
hulls <- mdsD %>%
  group_by(!!config$cov_factor) %>%
  group_modify(~ find_hull(.x)) %>%
  ungroup()

ph <- ggplot(mdsD,aes(x=dim1,y=dim2))
ph <- ph + geom_point(aes(colour=.data[[config$grouping_variable]]),size=4)
ph <- ph + geom_polygon(data=hulls, alpha=.2, aes(fill = .data[[config$cov_factor]]))
# ph <- ph + scale_fill_discrete("Status",labels=c("euthyroid","hyperthyroid"))
ph <- ph + labs(colour = "treatment")
ph <- ph + theme(
  legend.title=element_text(size=12),
  legend.text=element_text(size=10),
  axis.title=element_text(size=12),
  axis.text=element_text(size=10))
# ph

## UP TO HERE
## try with timepoint and weight as continuous, fecal score as categorical (or viceversa)

### CCA
writeLines(" - Canonical Correspondence Analysis")
H <- select(mdsD, all_of(continuous_covs))
preProcValues <- preProcess(H, method = c("knnImpute"))
H <- predict(preProcValues, H)

simX <- as.matrix(M)

vec <- sample(ncol(simX), nrow(M)+60)
X <- simX[,vec]

## safety net (filter to make sure all row sums are larger than zero)
vec <- rowSums(X) > 0
X <- X[vec,]
H <- H[vec,]

cca1 <- cca(X ~ timepoint+weight, data=H, scale=TRUE)

V <- as.data.frame(cca1$CCA$wa)
V$id <- row.names(V)
V <- inner_join(V,mdsD,by="id")

B <- as.data.frame(cca1$CCA$biplot)

pCCA <- ggplot(V,aes(CCA1,CCA2)) + geom_point(aes(colour = .data[[config$grouping_variable]],size=3))
pCCA <- pCCA + geom_segment(data = B, aes(xend = B[ ,"CCA1"], yend=B[ ,"CCA2"]),
                            x=0, y=0, colour="black",
                            arrow=arrow(angle=25, length=unit(0.5, "cm")))
pCCA <- pCCA + geom_text(data=B, aes(x=B[ ,"CCA1"], y=B[ ,"CCA2"], label=row.names(B)), 
                         size=6, vjust=1.5, colour="black")
pCCA <- pCCA + theme(strip.text = element_text(size=20),
                     axis.title=element_text(size=12),
                     axis.text=element_text(size=8))
pCCA <- pCCA + theme(legend.position="none")
# pCCA
# 
# multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#   require(grid)
#   
#   # Make a list from the ... arguments and plotlist
#   plots <- c(list(...), plotlist)
#   
#   numPlots = length(plots)
#   
#   # If layout is NULL, then use 'cols' to determine layout
#   if (is.null(layout)) {
#     # Make the panel
#     # ncol: Number of columns of plots
#     # nrow: Number of rows needed, calculated from # of cols
#     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#                      ncol = cols, nrow = ceiling(numPlots/cols))
#   }
#   
#   if (numPlots==1) {
#     print(plots[[1]])
#     
#   } else {
#     # Set up the page
#     grid.newpage()
#     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
#     
#     # Make each plot, in the correct location
#     for (i in 1:numPlots) {
#       # Get the i,j matrix positions of the regions that contain this subplot
#       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
#       
#       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
#                                       layout.pos.col = matchidx$col))
#     }
#   }
# }
# 
#  
# multiplot(pCCA,ph,cols=2)

writeLines(" - saving results")
fname = file.path(outdir, "cca_plot.png")
g <- ggarrange(pCCA, ph, nrow = 1, common.legend = FALSE, labels = c("A", "B"), widths = c(0.45,0.55))
ggsave(filename = fname, plot = g, device = "png", width = 9, height = 6, units = "in", dpi = 150)

fname = file.path(outdir, "cca_results.txt")
fwrite(x = V, file = fname)

print("DONE!!")

