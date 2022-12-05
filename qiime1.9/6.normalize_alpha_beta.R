## Rscript
## script to normalize the filtered OTU table

#############################################################################
## This script is mainly meant to be run locally
## where R packages can more easily be installed/updated
#############################################################################

## SET UP
library("ape")
library("ggplot2")
library("phyloseq")
library("tidyverse")
library("data.table")
library("metagenomeSeq")

## PARAMETERS
HOME <- Sys.getenv("HOME")
prj_folder = file.path(HOME, "Documents/USEFUL")
analysis_folder = "Analysis/microbiome_predipping/qiime1.9"
fname = "5.filter_OTUs/otu_table_filtered.biom"
conf_file = "Config/mapping_file.csv"
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
names(metadata)[1] <- "sample-id"
if(is.numeric(metadata$`sample-id`)) metadata$`sample-id` = paste("sample",metadata$`sample-id`,sep="-") # in case your sample-id are not only numeric, remove or comment if(is.numeric(metadata$`sample-id`))
metadata <- as.data.frame(metadata)
row.names(metadata) <- metadata$`sample-id`
metadata$`sample-id` <- NULL

## read into phyloseq
writeLines(" - add metadata to the phyloseq object")
samples = sample_data(metadata)
otu_tax_sample = phyloseq(otu,taxa,samples)
sample_data(otu_tax_sample)

# sample_names(otu_tax_sample)
# rank_names(otu_tax_sample)
# taxa_names(otu_tax_sample)

# plot_bar(otu_tax_sample, "Class", fill="Kingdom") + facet_grid(timepoint~treatment)

## making results folder
if(!file.exists(file.path(prj_folder, analysis_folder, "results"))) dir.create(file.path(prj_folder, analysis_folder, "results"), showWarnings = FALSE)

## making figures folder
if(!file.exists(file.path(prj_folder, analysis_folder, "results", "figures"))) dir.create(file.path(prj_folder, analysis_folder, "results", "figures"), showWarnings = FALSE)

## Alpha diversity
## alpha diversity is calculated on the original count data, not normalised 
## (see https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis)
## (see https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531)
writeLines(" - calculate alpha diversity indices")
alpha = estimate_richness(otu_tax_sample, split = TRUE)
alpha$"sample-id" = row.names(alpha)
alpha = relocate(alpha, `sample-id`)
fwrite(x = alpha, file = file.path(prj_folder, analysis_folder, "results", "alpha.csv"))
p <- plot_richness(otu_tax_sample, x="timepoint", color="treatment", shape = "herd")
ggsave(filename = file.path(prj_folder, analysis_folder, "results","figures", "alpha_plot.png"), plot = p, device = "png", width = 11, height = 7)

## Preprocessing: e.g. filtering
## making normalization folder
if(!file.exists(file.path(prj_folder, outdir))) dir.create(file.path(prj_folder, outdir), showWarnings = FALSE)

writeLines(" - CSS normalization")
otu_tax_sample_norm = phyloseq_transform_css(otu_tax_sample, norm = TRUE, log = FALSE)

## add taxonomy to normalised counts
otu_css_norm = base::as.data.frame(otu_table(otu_tax_sample_norm))
otu_css_norm$tax_id = row.names(otu_css_norm)
otu_css_norm <- relocate(otu_css_norm, tax_id)
taxonomy = as.data.frame(tax_table(otu_tax_sample_norm))
taxonomy$tax_id = row.names(taxonomy)
taxonomy <- relocate(taxonomy, tax_id)
otu_css_norm = otu_css_norm %>% inner_join(taxonomy, by = "tax_id")

writeLines(" - writing out the CSS normalized OTU table")
fwrite(x = otu_css_norm, file = file.path(prj_folder, outdir, "results", "otu_norm_CSS.csv"))

## relative abundances
otu_relative = transform_sample_counts(otu_tax_sample, function(x) x/sum(x) )
otu_rel_filtered = filter_taxa(otu_relative, function(x) mean(x) > 5e-3, TRUE)
nrow(otu_table(otu_rel_filtered))

writeLines(" - additionas plots")
random_tree = rtree(ntaxa((otu_rel_filtered)), rooted=TRUE, tip.label=taxa_names(otu_rel_filtered))
plot(random_tree)

biom1 = merge_phyloseq(otu_rel_filtered, random_tree)
plot_tree(biom1, color="treatment", label.tips="taxa_names", ladderize="left", plot.margin=0.3)
plot_tree(biom1, color="Genus", shape="timepoint", size="abundance")


###############
## distances ##
###############
writeLines(" - beta diversity: distance matrices")
writeLines(" - available distance metrics")
dist_methods <- unlist(distanceMethodList)
print(dist_methods)

## bray-curtis
writeLines(" - calculate Bray-Curtis distances")
distances = distance(otu_tax_sample_norm, method="bray", type = "samples")
iMDS  <- ordinate(otu_tax_sample_norm, "MDS", distance=distances)
p <- plot_ordination(otu_tax_sample_norm, iMDS, color="timepoint", shape="treatment")
ggsave(filename = file.path(prj_folder, analysis_folder, "results", "figures","mds_plot_bray_curtis.png"), plot = p, device = "png")

writeLines(" - write out distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, analysis_folder, "results", "bray_curtis_distances.csv"))

## euclidean
writeLines(" - calculate Euclidean distances")
distances = distance(otu_tax_sample_norm, method="euclidean", type = "samples")
iMDS  <- ordinate(otu_tax_sample_norm, "MDS", distance=distances)
p <- plot_ordination(otu_tax_sample_norm, iMDS, color="timepoint", shape="treatment")
ggsave(filename = file.path(prj_folder, analysis_folder, "results","figures", "mds_plot_euclidean.png"), plot = p, device = "png")

writeLines(" - write out euclidean distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, analysis_folder, "results", "euclidean_distances.csv"))

### add tree ###
random_tree = rtree(ntaxa((otu_tax_sample_norm)), rooted=TRUE, tip.label=taxa_names(otu_tax_sample_norm))
otu_norm_tree = merge_phyloseq(otu_tax_sample_norm, random_tree)
# plot_tree(otu_norm_tree, color="treatment", label.tips="taxa_names", ladderize="left", plot.margin=0.3)
plot_tree(otu_norm_tree, color="Phylum", shape="treatment", size="abundance")

## unifrac
writeLines(" - calculate Unifrac distances")
distances = distance(otu_norm_tree, method="unifrac", type = "samples")
iMDS  <- ordinate(otu_norm_tree, "MDS", distance=distances)
p <- plot_ordination(biom1, iMDS, color="timepoint", shape="treatment")
ggsave(filename = file.path(prj_folder, analysis_folder, "results", "figures","mds_plot_unifrac.png"), plot = p, device = "png")

writeLines(" - write out Unifrac distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, analysis_folder, "results", "unifrac_distances.csv"))

## weighted unifrac
writeLines(" - calculate weighted Unifrac distances")
distances = distance(otu_norm_tree, method="wunifrac", type = "samples")
iMDS  <- ordinate(otu_norm_tree, "MDS", distance=distances)
p <- plot_ordination(biom1, iMDS, color="timepoint", shape="treatment")
ggsave(filename = file.path(prj_folder, analysis_folder, "results", "figures","mds_plot_weighted_unifrac.png"), plot = p, device = "png")

writeLines(" - write out weighted Unifrac distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, analysis_folder, "results", "weighted_unifrac_distances.csv"))

