## Rscript
## script to normalize the filtered OTU table

#############################################################################
## This script is mainly meant to be run locally
## where R packages can more easily be installed/updated
#############################################################################

## SET UP
library("ape")
library("phyloseq")
library("tidyverse")
library("data.table")
library("metagenomeSeq")

## PARAMETERS
HOME <- Sys.getenv("HOME")
prj_folder = file.path(HOME, "Documents/USEFUL")
analysis_folder = "Analysis/microbiome_chlorine/micca"
fname = "otu_table/otu_table_filtered.biom"
conf_file = "Config/mapping_file_chlorine.csv"
min_tot_counts = 100 ## minimum number of total counts per sample to be included in the analysis
outdir = file.path(analysis_folder)

repo = file.path(prj_folder, "metabarcoding")
# source(file.path(prj_folder, repo, "r_scripts/dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
# source(file.path(prj_folder, repo, "r_scripts/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(repo, "support_functions/dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(repo, "support_functions/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/

writeLines(" - reading the filtered (OTU-wise) biom file into phyloseq")
## both the OTU table and the taxonomic classification are available from the biom file (qiime 1.9)
biom_otu_tax <- phyloseq::import_biom(BIOMfilename = file.path(prj_folder,analysis_folder,fname))

writeLines(" - removing samples with too few total counts")
biom_otu_tax = prune_samples(sample_sums(biom_otu_tax)>=min_tot_counts, biom_otu_tax)

otu = otu_table(biom_otu_tax, taxa_are_rows = TRUE)
taxa = tax_table(biom_otu_tax)

print(paste("N. of OTUs read from biom file is:", nrow(otu)))
print(paste("N .of samples retained after filtering is:", ncol(otu)))

colnames(otu) <- paste("sample-",colnames(otu),sep="")
print(head(otu))

writeLines(" - change the names of taxonomic levels to Kngdom, Class etc.")
colnames(taxa) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species") #if number does not fit, add "" as blank spaces to solve the problem
print(head(taxa))

## metadata
writeLines(" - reading the metadata")
metadata = fread(file.path(prj_folder,conf_file))
metadata = relocate(metadata, `sample-id`)
names(metadata)[1] <- "sample-id"
metadata$`sample-id` = paste("sample",metadata$`sample-id`,sep="-") 
if(is.numeric(metadata$`sample-id`)) metadata$`sample-id` = paste("sample",metadata$`sample-id`,sep="-") # in case your sample-id are not only numeric, remove or comment if(is.numeric(metadata$`sample-id`))
metadata <- as.data.frame(metadata)
row.names(metadata) <- metadata$`sample-id`
metadata$`sample-id` <- NULL

## read into phyloseq
writeLines(" - add metadata to the phyloseq object")
samples = sample_data(metadata)
otu_tax_sample = phyloseq(otu,taxa,samples)
sample_data(otu_tax_sample)

## save phyloseq object
save(otu_tax_sample, file = "Analysis/microbiome_chlorine/micca/otu_table/phyloseq.RData")
