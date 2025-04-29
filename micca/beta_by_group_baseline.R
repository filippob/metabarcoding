
## SET UP
library("ape")
library("knitr")
library("vegan")
library("ggpubr")
library("ggplot2")
library("phyloseq")
library("ggfortify")
library("tidyverse")
# library("speedyseq")
library("data.table")
library("metagenomeSeq")

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
    repo = "Documents/cremonesi/metabarcoding",
    prjfolder = "Documents/cremonesi/tamponi_vaginali",
    analysis_folder = "Analysis",
    conf_file = "Config/mapping_file.csv",
    suffix = "vaginal_swabs",
    nfactors = 1, ## n. of design variables (e.g. treatment and timpoint --> nfactors = 2)
    min_tot_n = 15,
    min_sample = 3,
    project = "",
    treatment_column = "treatment",
    subject_column = "CAMPIONI",
    force_overwrite = FALSE
  ))
}

HOME <- Sys.getenv("HOME")
repo = file.path(HOME, config$repo)
prjfolder = file.path(HOME, config$prjfolder)
outdir = file.path(prjfolder,config$analysis_folder)

fname = paste("beta_by_group_bsl-chg.config_",config$suffix,".RData", sep="")
fname = file.path(outdir, fname)
save(config, file = fname)

source(file.path(repo, "support_functions/dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(repo, "support_functions/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/

## loading data previously imported in phyloseq
fname = file.path(outdir, "phyloseq.RData")
load(fname)

## filter
min_tot_n = config$min_tot_n
min_sample = config$min_sample

# temp = filter_taxa(otu_tax_sample, function(x) sum(x > 1) > (0.20*length(x)), TRUE)
filtered_taxa = filter_taxa(otu_tax_sample, function(x) sum(x > 1) > min_sample, TRUE)
filtered_taxa = filter_taxa(filtered_taxa, function(x) sum(x) > min_tot_n, TRUE)
summary(rowSums(otu_table(filtered_taxa))) |> print()

## Preprocessing: e.g. filtering
## making normalization folder
if(!file.exists(file.path(outdir))) dir.create(file.path(outdir), showWarnings = FALSE)

writeLines(" - CSS normalization")
otu_tax_sample_norm = phyloseq_transform_css(filtered_taxa, norm = TRUE, log = FALSE)
sample_data(otu_tax_sample_norm) |> head() |> print()

## subset for subproject
if(nchar(config$project) > 0) {
  
  otu_norm_subset <- subset_samples(otu_tax_sample_norm, project == config$project)
} else otu_norm_subset <- otu_tax_sample_norm

## changing treatment column
x <- sample_data(otu_tax_sample_norm)
sample_data(otu_norm_subset)$treatment <- dplyr::pull(x, !!config$treatment_column )
sample_data(otu_norm_subset) |> head() |> print()

## function to convert otu table from phyloseq object to (transposed) matrix
vegan_otu <-  function(physeq){
  OTU <-  phyloseq::otu_table(physeq)
  if(phyloseq::taxa_are_rows(OTU)){
    OTU <-  t(OTU)
  }
  return(as(OTU,"matrix"))
}

steps = unique(sample_data(otu_norm_subset)$timepoint)

res = list(NULL)

for (k in steps) {
  
  print(k)
  
  otu <- subset_samples(otu_norm_subset, timepoint == k)
  
  temp <- vegan_otu(otu)
  ids = rownames(temp)
  temp <- as_tibble(temp)
  temp$sample_id = ids
  
  metadata = sampleidsmetadata = sample_data(otu)
  metadata$sample_id = row.names(metadata)
  
  temp <- temp |> inner_join(metadata, by = "sample_id")
  temp <- temp |> relocate(c("sample_id", "CAMPIONI","timepoint","treatment"), .before = "DENOVO1")
  temp <- temp |> arrange(CAMPIONI)
  
  res[[k]] = temp
  
}

otu_chg <- (as.matrix(res$T2[,-c(1:4)]) - as.matrix(res$T1[,-c(1:4)]))
otu_chg <- cbind(res$T2[,c(1:4)], otu_chg)

###############
## distances ##
###############
writeLines(" - beta diversity: distance matrices")
writeLines(" - available distance metrics")

# Calculate Bray-Curtis dissimilarity
x <- as.matrix(otu_chg[,-c(1:4)])

## normalise matrix in 0-1
x <- 100 * (x-min(x))/(max(x)-min(x)) 
## normalise matrix in 0-1
# x <- (x-mean(x))/(sd(x))
# x <- x + abs(min(x))

bray_curtis_dist <- vegdist(x, method = "bray", binary = FALSE)

# Creating easy to view matrix and writing .csv
distmat <- as.matrix(bray_curtis_dist, labels = TRUE)
similarity_matrix = 1-distmat
rownames(similarity_matrix) <- paste(otu_chg$treatment, otu_chg$CAMPIONI, sep="_")
colnames(similarity_matrix) <- paste(otu_chg$treatment, otu_chg$CAMPIONI, sep="_")

fname = paste("distances_bsl_chg", config$suffix, ".png", sep="")
fname = file.path(outdir, "figures", fname)

png(fname)
heatmap(similarity_matrix)
dev.off()

fpath = file.path(outdir, "tables")
dir.create(fpath, showWarnings = FALSE)
fname = file.path(fpath, paste("bray_curtis_",config$suffix,".csv",sep=""))
fwrite(x = distmat, file = fname)

### multidimensional scaling ##
mds_D <- distmat %>%
  cmdscale(k=3) %>%
  as_tibble() %>%
  rename(dim1 = V1, dim2 = V2, dim3 = V3) %>%
  mutate(id = otu_chg$CAMPIONI, treatment = otu_chg$treatment)

g <- ggscatter(data = mds_D, x = "dim1", y = "dim2",
               label = NULL,
               color = "treatment",
               palette = "jco",
               size = 1, 
               ellipse = TRUE,
               ellipse.type = "convex",
               repel = TRUE)

g

fname = paste("mds_plot_bray-curtis_ellipse_bsl-chg_", config$suffix, "_", "treatment" , ".png")
to_save[[paste("beta_div_plot", "treatment", sep="_")]] = g
ggsave(filename = file.path(prjfolder, config$analysis_folder, "figures", fname), plot = g, device = "png")


## 3D plot
library("plotly")

p <- plot_ly(data = mds_D, 
             x = ~dim1, y = ~dim2, z = ~dim3,
             type = "scatter3d",
             color = ~treatment,
             size = 1) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC 1'),
                      yaxis = list(title = 'PC 2'),
                      zaxis = list(title = 'PC 3')),
         annotations = list(
           x = 0.005,
           y = 0.01,
           text = 'Vaginal Swabs',
           xref = 'pc1',
           yref = 'pc2',
           showarrow = FALSE
         ))

# p

#####################
## PERMANOVA
#####################

print("Permanova on the whole dataset: comparison between types, treatmetns and type x treatment:")
source("~/Documents/cremonesi/rumine_anafi/parwise.adonis.r")

fpath = file.path(outdir, "tables")
dir.create(fpath, showWarnings = FALSE)
fname = file.path(fpath, paste("permanova_bsl-chg",config$suffix,".csv",sep=""))
fwrite(x = list("PERMANOVA"), file = fname)

write(x = "PERMANOVA", file = fname)

temp <- cbind(distmat, otu_chg[,4])

if (config$nfactors > 1) {
  
  nvars = config$nfactors
  matx= data.matrix(temp[,seq(1,ncol(temp)-nvars)])
  
  obj <- adonis2(matx ~ otu_chg$timepoint+otu_chg$treatment, permutations = 1000)
  print(kable(obj))
  # fwrite(x = list(kable(obj)), file = fname, append = TRUE)
  write(x = "ALL DATA", file = fname, append = TRUE)
  write(x = kable(obj), file = fname, append = TRUE)
  
  to_save <- list("permanova_all"=obj)
  
  tipos = unique(otu_chg$timepoint)
  for (tt in tipos) {
    
    print(tt)
    vec = (otu_chg$timepoint == tt)
    temp2 <- otu_chg[vec,c(TRUE,vec,TRUE,TRUE)]
    temp1 <- matx[vec,vec]
    
    print(paste("Permanova on the", tt, "subset: comparison between treatments:"))
    obj <- adonis2(temp1 ~ temp2$treatment, permutations = 1000)
    write(x = tt, file = fname, append = TRUE)
    write(x = kable(obj), file = fname, append = TRUE)
    print(kable(obj))
    
    to_save[[paste("permanova", tt, sep="_")]] = obj
  }
} else {
  
  nvars = config$nfactors
  matx= data.matrix(temp[,seq(1,ncol(temp)-nvars)])
  
  obj <- adonis2(matx ~ otu_chg$treatment, permutations = 1000)
  print(kable(obj))
  write(x = kable(obj), file = fname, append = TRUE)
  
  to_save <- list("permanova"=obj)
}


