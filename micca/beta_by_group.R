
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
    prjfolder = "Documents/Suikerbiet/its_2025",
    analysis_folder = "Analysis",
    conf_file = "Config/mapping_file.csv",
    suffix = "its1f-its2",
    nfactors = 2, ## n. of design variables (e.g. treatment and timpoint --> nfactors = 2)
    min_tot_n = 15,
    min_sample = 3,
    covariates = "", ## string with covariates separated by a comma
    project = "",
    treatment_column = "treatment",
    timepoint_column = "timepoint",
    force_overwrite = FALSE
  ))
}

HOME <- Sys.getenv("HOME")
repo = file.path(HOME, config$repo)
prjfolder = file.path(HOME, config$prjfolder)
outdir = file.path(prjfolder,config$analysis_folder)

fname = paste("beta_by_group.config_",config$suffix,".RData", sep="")
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

steps = sample_data(otu_norm_subset) |>
  pull(!!config$timepoint_column) |>
  unique()

###############
## distances ##
###############
writeLines(" - beta diversity: distance matrices")
writeLines(" - available distance metrics")

if (all(steps == "") || is.null(steps)) {
  
  print(paste(" - calculate Bray-Curtis distances for ", config$treatment_column))
  distances = distance(otu_norm_subset, method="bray", type = "samples")
  iMDS  <- ordinate(otu_norm_subset, "MDS", distance=distances)
  p <- plot_ordination(otu_norm_subset, iMDS, color="treatment")
  fname = paste("mds_plot_bray_curtis_", config$suffix, "_treatment.png")
  ggsave(filename = file.path(prjfolder, config$analysis_folder, "figures", fname), plot = p, device = "png")
}

## bray-curtis
for (k in steps) {
  
  print(paste(" - calculate Bray-Curtis distances for ", k))
  temp = subset_samples(otu_norm_subset, timepoint == k)
  distances = distance(temp, method="bray", type = "samples")
  iMDS  <- ordinate(temp, "MDS", distance=distances)
  p <- plot_ordination(temp, iMDS, color="treatment")
  fname = paste("mds_plot_bray_curtis_", config$suffix, k , ".png")
  ggsave(filename = file.path(prjfolder, config$analysis_folder, "figures", fname), plot = p, device = "png")
}

## permanova
metadata <- sample_data(otu_norm_subset)
metadata$`sample-id` = row.names(metadata)
row.names(metadata) <- NULL
metadata <- as_tibble(metadata)

if (config$covariates != "") { covariates = unlist(strsplit(config$covariates, split = ",")) 
} else covariates = NULL
print(paste("The following covariates are used:", covariates))

if (config$treatment_column != "treatment") metadata$treatment = NULL

metadata <- metadata |> rename('sample-id' = !!config$sample_column, treatment = !!config$treatment_column)

if("timepoint" %in% names(metadata)) {
  metadata = metadata |> select(c(`sample-id`, treatment, timepoint, all_of(covariates))) |> mutate(treatment = as.factor(treatment))
} else metadata = metadata |> select(c(`sample-id`, treatment, all_of(covariates))) |> mutate(treatment = as.factor(treatment))

## UNCOMMENT BELOW IF YOU NEED TO CHANGE TREATMENT LABELS
# old_treat = unique(sample_data(otu_norm_subset)$treatment)
# new_treat = c("non-EU-", "PC", "EU", "non-EU+")
# sample_data(otu_norm_subset)$treatment = new_treat[match(sample_data(otu_norm_subset)$treatment, old_treat)]
# sample_data(otu_norm_subset)$treatment

distances = distance(otu_norm_subset, method="bray", type = "samples")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")

# temp = dplyr::select(metadata, c(`sample-id`,timepoint,Antibiotic)) |> rename(treatment = Antibiotic)
if(config$nfactors > 1) {
  
  temp = dplyr::select(metadata, c(`sample-id`, timepoint, treatment, all_of(covariates)))
} else temp = dplyr::select(metadata, c(`sample-id`,treatment, all_of(covariates)))

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
  
  nvars = config$nfactors + length(covariates)
  matx = data.matrix(temp[,seq(1,ncol(temp)-nvars)])
  
  if(config$covariates != "") { 
    covar = paste(gsub(",", "+", config$covariates))
    obj <- adonis2(matx ~ dx$timepoint+dx$treatment + dx[[covar]], permutations = 1000)
    } else obj <- adonis2(matx ~ dx$timepoint+dx$treatment, permutations = 1000)
  
  print(kable(obj))
  # fwrite(x = list(kable(obj)), file = fname, append = TRUE)
  write(x = "ALL DATA", file = fname, append = TRUE)
  write(x = kable(obj), file = fname, append = TRUE)
  
  to_save <- list("permanova_all"=obj)
  
  tipos = unique(dx$timepoint)
  for (tt in tipos) {
    
    print(tt)
    vec = (dx$timepoint == tt)
    temp2 <- dx[vec,c(TRUE,vec,rep(TRUE,nvars))]
    temp1 <- matx[vec,vec]
    
    print(paste("Permanova on the", tt, "subset: comparison between treatments:"))
    if(config$covariates != "") { 
      obj <- adonis2(temp1 ~ temp2$treatment + temp2[[covar]], permutations = 1000)
    } else obj <- adonis2(temp1 ~ temp2$treatment, permutations = 1000)
    
    write(x = tt, file = fname, append = TRUE)
    write(x = kable(obj), file = fname, append = TRUE)
    print(kable(obj))
    
    to_save[[paste("permanova", tt, sep="_")]] = obj
  }
} else {
  
  covar = paste(gsub(",", "+", config$covariates))
  nvars = config$nfactors
  matx= data.matrix(temp[,seq(1,ncol(temp)-nvars)])
  
  obj <- adonis2(matx ~ dx$treatment + dx[[covar]], permutations = 1000)
  print(kable(obj))
  write(x = kable(obj), file = fname, append = TRUE)
  
  to_save <- list("permanova"=obj)
}


library("ggpubr")
library("ggfortify")

print(paste("processing stratifying variable",config$treatment_column))
mtd <- sample_data(otu_norm_subset)
distances = distance(otu_norm_subset, method="bray", type = "samples")
iMDS  <- ordinate(otu_norm_subset, "MDS", distance=distances)
mds_subset <- as_tibble(iMDS$vectors)
mds_subset$treatment = mtd$treatment
mds_subset$timepoint = mtd$timepoint

gall <- ggscatter(data = mds_subset, x = "Axis.1", y = "Axis.2",
               label = NULL,
               color = "timepoint",
               shape = "treatment",
               palette = "jco",
               size = 3,
               ellipse = TRUE,
               ellipse.type = "norm",
               repel = TRUE) + scale_shape_manual(values = c(4,19,2)) 

fname = paste("mds_plot_bray-curtis_ellipse_", config$suffix, "_", "ALL" , ".png")
to_save[[paste("beta_div_plot", "treatment", sep="_")]] = gall
ggsave(filename = file.path(prjfolder, config$analysis_folder, "figures", fname), plot = gall, device = "png", dpi = 150, height = 6, width = 7)



if (all(steps == "") || is.null(steps)) {
  
  print(paste("processing stratifying variable",config$treatment_column))
  mtd <- sample_data(otu_norm_subset)
  distances = distance(otu_norm_subset, method="bray", type = "samples")
  iMDS  <- ordinate(otu_norm_subset, "MDS", distance=distances)
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
  
  fname = paste("mds_plot_bray-curtis_ellipse_", config$suffix, "_", "treatment" , ".png")
  to_save[[paste("beta_div_plot", "treatment", sep="_")]] = g
  ggsave(filename = file.path(prjfolder, config$analysis_folder, "figures", fname), plot = g, device = "png")
}

for (k in steps) {
  
  print(paste("processing stratifying variable",k))
  temp <- subset_samples(otu_norm_subset, timepoint == k)
  mtd <- sample_data(temp)
  mtd$treatment = c("0","1","2")[match(mtd$treatment,unique(mtd$treatment))]
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
                 repel = TRUE) + theme(legend.title=element_blank())
  
  fname = paste("mds_plot_bray-curtis_ellipse_", config$suffix, "_", k , ".png")
  to_save[[paste("beta_div_plot", k, sep="_")]] = g
  ggsave(filename = file.path(prjfolder, config$analysis_folder, "figures", fname), plot = g, device = "png")
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
                 repel = TRUE) + theme(legend.title=element_blank())
  
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

