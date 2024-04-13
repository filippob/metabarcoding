
## SET UP
library("ape")
library("knitr")
library("vegan")
library("ggplot2")
library("phyloseq")
library("tidyverse")
# library("speedyseq")
library("data.table")
library("metagenomeSeq")

## PARAMETERS
args = commandArgs(trailingOnly=TRUE)
if (length(args) >= 1) {
  
  #loading the parameters
  source(args[1])
  # source("Analysis/hrr/config.R")
  
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
    prjfolder = "Documents/moroni/capre/delower",
    analysis_folder = "Analysis/results",
    conf_file = "Config/mapping_file.csv",
    suffix = "goat_milk",
    nfactors = 2, ## n. of design variables (e.g. treatment and timpoint --> nfactors = 2)
    min_tot_n = 15,
    min_sample = 3,
    project = "",
    treatment_column = "Antibiotic",
    force_overwrite = FALSE
  ))
}

HOME <- Sys.getenv("HOME")
repo = file.path(HOME, config$repo)
prjfolder = file.path(HOME, config$prjfolder)
outdir = file.path(prjfolder,config$analysis_folder)

fname = file.path(outdir, "beta_by_group.config.r")
fwrite(x = config, file = fname)

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
sample_data(otu_norm_subset)$treatment <- dplyr::pull(x, !!config$treatment_column )
sample_data(otu_norm_subset) |> head() |> print()

steps = unique(sample_data(otu_norm_subset)$timepoint)

###############
## distances ##
###############
writeLines(" - beta diversity: distance matrices")
writeLines(" - available distance metrics")
## bray-curtis
for (k in steps) {
  
  print(paste(" - calculate Bray-Curtis distances for ", k))
  temp = subset_samples(otu_norm_subset, timepoint == k)
  distances = distance(temp, method="bray", type = "samples")
  iMDS  <- ordinate(temp, "MDS", distance=distances)
  p <- plot_ordination(temp, iMDS, color="treatment")
  fname = paste("mds_plot_bray_curtis_", config$suffix, k , ".png")
  ggsave(filename = file.path(prj_folder, analysis_folder, "figures", fname), plot = p, device = "png")
}

## permanova
metadata <- sample_data(otu_norm_subset)
metadata$`sample-id` = row.names(metadata)
row.names(metadata) <- NULL
metadata <- as_tibble(metadata)

distances = distance(otu_norm_subset, method="bray", type = "samples")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")

# temp = dplyr::select(metadata, c(`sample-id`,timepoint,Antibiotic)) |> rename(treatment = Antibiotic)
temp = dplyr::select(metadata, c(`sample-id`,timepoint,treatment))
dx <- dx %>% inner_join(temp, by = c("row" = "sample-id"))

temp = select(dx,-row)

print("Permanova on the whole dataset: comparison between types, treatmetns and type x treatment:")
source("~/Documents/cremonesi/rumine_anafi/parwise.adonis.r")

fpath = file.path(outdir, "tables")
dir.create(fpath, showWarnings = FALSE)
fname = file.path(fpath, paste("permanova_",config$suffix,".csv",sep=""))
fwrite(x = list("PERMANOVA"), file = fname)

write(x = "PERMANOVA", file = fname)

if (config$nfactors > 1) {
  
  nvars = config$nfactors
  matx= data.matrix(temp[,seq(1,ncol(temp)-nvars)])
  
  obj <- adonis2(matx ~ dx$timepoint+dx$treatment, permutations = 1000)
  print(kable(obj))
  # fwrite(x = list(kable(obj)), file = fname, append = TRUE)
  write(x = "ALL DATA", file = fname, append = TRUE)
  write(x = kable(obj), file = fname, append = TRUE)
  
  to_save <- list("permanova_all"=obj)
  
  tipos = unique(dx$timepoint)
  for (tt in tipos) {
    
    print(tt)
    vec = (dx$timepoint == tt)
    temp2 <- dx[vec,c(TRUE,vec,TRUE,TRUE)]
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
  
  obj <- adonis2(matx ~ dx$treatment, permutations = 1000)
  print(kable(obj))
  write(x = kable(obj), file = fname, append = TRUE)
  
  to_save <- list("permanova"=obj)
}


library("ggpubr")
library("ggfortify")

for (k in steps) {
  
  print(paste("processing stratifying variable",k))
  temp <- subset_samples(otu_norm_subset, timepoint == k)
  mtd <- sample_data(temp)
  distances = distance(temp, method="bray", type = "samples")
  iMDS  <- ordinate(temp, "MDS", distance=distances)
  mds_subset <- as_tibble(iMDS$vectors)
  mds_subset$treatment = mtd$treatment
  
  g <- ggscatter(data = mds_subset, x = "Axis.1", y = "Axis.2",
                 label = NULL,
                 color = "treatment",
                 palette = "jco",
                 size = 1,
                 ellipse = TRUE,
                 ellipse.type = "norm",
                 repel = TRUE)
  
  fname = paste("mds_plot_bray-curtis_ellipse_", config$suffix, "_", k , ".png")
  to_save[[paste("beta_div_plot", k, sep="_")]] = g
  ggsave(filename = file.path(prj_folder, analysis_folder, "figures", fname), plot = g, device = "png")
}


if(length(steps) > 1) {
  
  mtd <- sample_data(otu_norm_subset)
  distances = distance(otu_norm_subset, method="bray", type = "samples")
  iMDS  <- ordinate(temp, "MDS", distance=distances)
  mds_subset <- as_tibble(iMDS$vectors)
  mds_subset$timepoint = mtd$timepoint
  g <- ggscatter(data = mds_subset, x = "Axis.1", y = "Axis.2",
                 label = NULL,
                 color = "timepoint",
                 palette = "jco",
                 size = 1,
                 ellipse = TRUE,
                 ellipse.type = "norm",
                 repel = TRUE)
  
  to_save[["beta_div_plot_timepoint"]] = g
  
  ## distance matrix (Bray-Curtis, all dataset)
  to_save[["distance_matrix"]] = distances
}

g

## save results to R object
fname = paste("beta_results_", config$suffix, ".RData", sep="")
fname = file.path(outdir, fname)
save(to_save, file = fname)

print("DONE!")
