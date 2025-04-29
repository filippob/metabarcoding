## Rscript
## script to normalize the filtered OTU table

#############################################################################
## This script is mainly meant to be run locally
## where R packages can more easily be installed/updated
#############################################################################

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("metagenomeSeq")

## SET UP
library("ape")
library("tools")
library("phyloseq")
library("tidyverse")
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
    #base_folder = '~/Documents/SMARTER/Analysis/hrr/',
    #genotypes = "Analysis/hrr/goat_thin.ped",
    repo = "Documents/cremonesi/metabarcoding",
    prjfolder = "Documents/moroni/lettiera",
    analysis_folder = "Analysis",
    fname = "filtered_otu/otu_table_filtered.biom",
    conf_file = "Config/mapping_file.csv",
    suffix = "bedding_wk2",
    nfactors = 2, ## n. of design variables (e.g. treatment and timpoint --> nfactors = 2)
    min_tot_n = 10,
    min_sample = 2,
    project = "", ## USE ONLY FOR SUBSETTING !!
    treatment_column = "treatment",
    sample_column = "sample_id",
    subset_variable = "week",
    subset_group = "WEEK2", ## subset data by sample variable (e.g. experiment, group, sex, etc.),
    force_overwrite = FALSE
  ))
}

HOME <- Sys.getenv("HOME")
repo = file.path(HOME, config$repo)
prjfolder = file.path(HOME, config$prjfolder)
analysis_folder = file.path(prjfolder,config$analysis_folder)
fname = config$fname

config_fname = file.path(analysis_folder, "import_phyloseq.config.RData")
save(config, file = config_fname)

# source(file.path(prj_folder, repo, "r_scripts/dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
# source(file.path(prj_folder, repo, "r_scripts/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(repo, "support_functions/dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(repo, "support_functions/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/

writeLines(" - reading the filtered (OTU-wise) biom file into phyloseq")
## both the OTU table and the taxonomic classification are available from the biom file (qiime 1.9)
biom_otu_tax <- phyloseq::import_biom(BIOMfilename = file.path(analysis_folder,fname))
# otu = otu_table(biom_otu_tax, taxa_are_rows = TRUE)
# ncol(otu)

writeLines(" - removing samples with too few total counts")
biom_otu_tax = prune_samples(sample_sums(biom_otu_tax) >= config$min_tot_n, biom_otu_tax)

otu = otu_table(biom_otu_tax, taxa_are_rows = TRUE)
taxa = tax_table(biom_otu_tax)

print(paste("N. of OTUs read from biom file is:", nrow(otu)))
print(paste("N .of samples retained after filtering is:", ncol(otu)))

colnames(otu) <- paste("sample-",colnames(otu),sep="")
print(head(otu))

writeLines(" - change the names of taxonomic levels to Kngdom, Class etc.")
colnames(taxa) <- c("Kingdom","Phylum","Class","Order","Family","Genus") #if number does not fit, add "" as blank spaces to solve the problem
# colnames(taxa) <- c("Kingdom","Phylum","Class","Order","Family","Genus") #from RDP
print(head(taxa))

## metadata
writeLines(" - reading the metadata")
metadata = fread(file.path(prjfolder, config$conf_file))
#metadata = metadata |> rename(`sample-id` = id) |> relocate(`sample-id`)
names(metadata)[1] <- "sample-id"
if (config$treatment_column != "treatment") metadata <- rename(metadata, 'treatment' = !!config$treatment_column)
metadata$`sample-id` = paste("sample",metadata$`sample-id`,sep="-") 
if(is.numeric(metadata$`sample-id`)) metadata$`sample-id` = paste("sample",metadata$`sample-id`,sep="-") # in case your sample-id are not only numeric, remove or comment if(is.numeric(metadata$`sample-id`))
metadata <- as.data.frame(metadata)
row.names(metadata) <- metadata$`sample-id`
metadata$`sample-id` <- NULL
if("timepoint" %in% names(metadata)) metadata$timepoint = as.factor(metadata$timepoint)
metadata$treatment = as.factor(metadata$treatment)

## read into phyloseq
writeLines(" - add metadata to the phyloseq object")
samples = sample_data(metadata)
otu_tax_sample = phyloseq(otu,taxa,samples)
sample_data(otu_tax_sample) |> head()
sample_data(otu_tax_sample) |> nrow()

## subset data if needed
if (!(is.null(config$treatment_column) | config$subset_variable == "") & !(is.null(config$subset_group) | config$subset_group == "")) {
  
  subset_group = config$subset_group
  print(paste("subsetting data by", subset_group))
  otu_tax_sample <- subset_samples(otu_tax_sample, get(eval(config$subset_variable)) == subset_group)
  # otu_tax_sample <- subset_samples(otu_tax_sample, week == subset_group)
  print("n. of samples left after subsetting")
  sample_data(otu_tax_sample) |> nrow() |> print()
}

## remove samples if treatment or timepoint is missing
if("timepoint" %in% names(metadata)) otu_tax_sample <- subset_samples(otu_tax_sample, !(is.na(timepoint)))
otu_tax_sample <- subset_samples(otu_tax_sample, !(is.na(treatment)))
sample_data(otu_tax_sample) |> nrow()

## save phyloseq object
dir.create(file.path(analysis_folder), showWarnings = FALSE)
fname = file.path(analysis_folder, "phyloseq.RData")
save(otu_tax_sample, file = fname)

