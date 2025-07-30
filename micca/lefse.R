## libraries
library("dplyr")
library("lefser")
library("ggplot2")
library("phyloseq")
library("SummarizedExperiment")

## PARAMETERS
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 1) {
  
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
    prjfolder = "Documents/Suikerbiet/its_2025/its_kyo",
    analysis_folder = "Analysis",
    conf_file = "Config/mapping_file.csv",
    suffix = "its-kyo",
    nfactors = 2, ## n. of design variables (e.g. treatment and timpoint --> nfactors = 2)
    min_tot_n = 15,
    min_sample = 3,
    project = "",
    sample_column = "sample-id",
    treatment_column = "treatment",
    factor_cov = "timepoint",
    sample_id = "", ## from metadata
    force_overwrite = FALSE
  ))
}

## SET UP
HOME <- Sys.getenv("HOME")
repo = file.path(HOME, config$repo)
prjfolder = file.path(HOME, config$prjfolder)
outdir = file.path(prjfolder,config$analysis_folder)

fname = paste("lefse.config_",config$suffix,".RData", sep="")
fname = file.path(outdir, fname)
save(config, file = fname)

source(file.path(HOME, config$repo, "support_functions/dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(HOME, config$repo, "support_functions/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMis

## loading data previously imported in phyloseq
fname = file.path(outdir, "phyloseq.RData")
load(fname)

otu_tax_sample_norm = phyloseq_transform_css(otu_tax_sample, norm = TRUE, log = FALSE)
otu_table = otu_tax_sample_norm

rm(otu_tax_sample)
rm(otu_tax_sample_norm)
gc()

counts <- unclass(otu_table(otu_table))
metadata <- as(sample_data(otu_table), "data.frame")
taxons <- unclass(tax_table(otu_table))

## PERFORM LEFSE PER STEP OF SECONDARY FACTOR VARIABLE (E.G. TIMEPOINT)
steps = pull(metadata, !!config$factor_cov) |> unique()
for (k in steps) {
  
  print(k)
  exp_metadata <- filter(metadata, .data[[config$factor_cov]] == k)
  exp_metadata$sample = rownames(exp_metadata)
  
  ids = row.names(exp_metadata)
  exp_df = counts |>
    as_tibble() |>
    select(all_of(ids)) |>
    as.matrix()
  
  rownames(exp_df) <- rownames(counts)
  print(paste("N. of OTUs is:", nrow(exp_df)))
  
  ## Create a SummarizedExperiment object
  exp_data = SummarizedExperiment(assays = list(counts = exp_df), colData = exp_metadata)
  tn <- get_terminal_nodes(rownames(exp_data))
  exp_data_n <- exp_data[tn,]
  exp_data_ra <- relativeAb(exp_data)
  
  print(" running LefSe model")
  res <- lefser(relab = exp_data_ra, classCol = "treatment")

  genus = taxons |>
    as_tibble() |>
    pull(Genus)

  otu_ids = row.names(taxons)
  res$features = genus[match(res$features,otu_ids)]
  print(head(res))
  
  fname = paste("lefse_", k, ".csv", sep="")
  path = file.path(outdir, "lefse")
  dir.create(path, showWarnings = FALSE)
  fname = file.path(path, fname)
  write.csv(x = na.omit(res), file = fname)
  
  fname = paste("lefse_", k, ".png", sep="")
  fname = file.path(path, fname)
  
  n = nrow(na.omit(res))
  
  print(paste("saving image to file", fname))
  plot_obj <- lefserPlot(na.omit(res))
  # class(plot_obj)
  ggsave(fname, plot = plot_obj, width = 6, height = n*0.5, dpi = 250, bg = "white")
}

print("DONE!!")
