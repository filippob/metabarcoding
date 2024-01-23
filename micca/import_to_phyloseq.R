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
library("phyloseq")
library("tidyverse")
library("data.table")
library("metagenomeSeq")

## PARAMETERS
HOME <- Sys.getenv("HOME")
repo = file.path(HOME, "Documents/cremonesi/metabarcoding")
prj_folder = file.path(HOME, "Documents/moroni/capre/delower")
analysis_folder = "Analysis"
fname = "filtered_otu/otu_table_filtered.biom"
conf_file = "Config/mapping_file.csv"
min_tot_counts = 15 ## minimum number of total counts per sample to be included in the analysis
outdir = file.path(analysis_folder,"results")
subset_group = "" ## subset data by sample variable (e.g. experiment, group, sex, etc.)

# source(file.path(prj_folder, repo, "r_scripts/dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
# source(file.path(prj_folder, repo, "r_scripts/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(repo, "support_functions/dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(repo, "support_functions/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/

writeLines(" - reading the filtered (OTU-wise) biom file into phyloseq")
## both the OTU table and the taxonomic classification are available from the biom file (qiime 1.9)
biom_otu_tax <- phyloseq::import_biom(BIOMfilename = file.path(prj_folder,analysis_folder,fname))

writeLines(" - removing samples with too few total counts")
biom_otu_tax = prune_samples(sample_sums(biom_otu_tax) >= min_tot_counts, biom_otu_tax)

otu = otu_table(biom_otu_tax, taxa_are_rows = TRUE)
taxa = tax_table(biom_otu_tax)

print(paste("N. of OTUs read from biom file is:", nrow(otu)))
print(paste("N .of samples retained after filtering is:", ncol(otu)))

colnames(otu) <- paste("sample-",colnames(otu),sep="")
print(head(otu))

writeLines(" - change the names of taxonomic levels to Kngdom, Class etc.")
colnames(taxa) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species") #if number does not fit, add "" as blank spaces to solve the problem
# colnames(taxa) <- c("Kingdom","Phylum","Class","Order","Family","Genus") #from RDP
print(head(taxa))

## metadata
writeLines(" - reading the metadata")
metadata = fread(file.path(prj_folder,conf_file))
#metadata = metadata |> rename(`sample-id` = id) |> relocate(`sample-id`)
names(metadata)[1] <- "sample-id"
metadata$`sample-id` = paste("sample",metadata$`sample-id`,sep="-") 
if(is.numeric(metadata$`sample-id`)) metadata$`sample-id` = paste("sample",metadata$`sample-id`,sep="-") # in case your sample-id are not only numeric, remove or comment if(is.numeric(metadata$`sample-id`))
metadata <- as.data.frame(metadata)
row.names(metadata) <- metadata$`sample-id`
metadata$`sample-id` <- NULL
metadata$timepoint = as.factor(metadata$timepoint)
metadata$treatment = as.factor(metadata$treatment)

## read into phyloseq
writeLines(" - add metadata to the phyloseq object")
samples = sample_data(metadata)
otu_tax_sample = phyloseq(otu,taxa,samples)
sample_data(otu_tax_sample) |> head()

## subset data if needed
if (!(is.null(subset_group) | subset_group == "")) {
  
  print(paste("subsetting data by", subset_group))
  otu_tax_sample <- subset_samples(otu_tax_sample, experiment == subset_group)
  print("n. of samples left after subsetting")
  sample_data(otu_tax_sample) |> nrow() |> print()
}

## remove samples if treatment or timepoint is missing
otu_tax_sample <- subset_samples(otu_tax_sample, !(is.na(treatment) | is.na(timepoint)))

## save phyloseq object
dir.create(file.path(prj_folder, outdir), showWarnings = FALSE)
fname = file.path(prj_folder, outdir, "phyloseq.RData")
save(otu_tax_sample, file = fname)

