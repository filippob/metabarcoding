## Rscript
## script to normalize the filtered OTU table

#############################################################################
## This script is mainly meant to be run locally
## where R packages can more easily be installed/updated
#############################################################################

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("phyloseq")

## SET UP
library("ape")
library("ggplot2")
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
    prjfolder = "Documents/Suikerbiet/its_2025",
    analysis_folder = "Analysis",
    conf_file = "Config/mapping_file.csv",
    suffix = "its1f-its2",
    nfactors = 2, ## n. of design variables (e.g. treatment and timpoint --> nfactors = 2)
    min_tot_n = 20,
    min_sample = 3,
    project = "",
    treatment_column = "treatment",
    sample_column = "sample",
    grouping_variable2 = "timepoint",
    grouping_variable1 = "treatment",
    force_overwrite = FALSE
  ))
}

HOME <- Sys.getenv("HOME")
repo = file.path(HOME, config$repo)
prjfolder = file.path(HOME, config$prjfolder)
outdir = file.path(prjfolder,config$analysis_folder)

config_fname = file.path(outdir, "normalise_diversity.config.RData")
save(config, file = config_fname)

## treatment levels as in the metadata file
grouping_variable1 = config$grouping_variable1
grouping_variable2 = config$grouping_variable2

repo = file.path(HOME, "Documents/cremonesi/metabarcoding")
source(file.path(repo, "support_functions/dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(repo, "support_functions/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/

## loading data previously imported in phyloseq
fname = file.path(outdir, "phyloseq.RData")
load(fname)

## making results folder
# if(!file.exists(file.path(prj_folder, analysis_folder, "results"))) dir.create(file.path(prj_folder, analysis_folder, "results"), showWarnings = FALSE)

## making figures folder
if(!file.exists(file.path(outdir, "figures"))) dir.create(file.path(outdir, "figures"), showWarnings = FALSE)

## Alpha diversity
## alpha diversity is calculated on the original count data, not normalised 
## (see https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis)
## (see https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531)
writeLines(" - calculate alpha diversity indices")
## alpha div estimated on not normalised OTU counts: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531
alpha = estimate_richness(otu_tax_sample, split = TRUE)
alpha$"sample-id" = row.names(alpha)
alpha = relocate(alpha, `sample-id`)
fwrite(x = alpha, file = file.path(outdir, "alpha.csv"))

## Preprocessing: e.g. filtering
## making normalization folder
if(!file.exists(file.path(outdir))) dir.create(file.path(outdir), showWarnings = FALSE)

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
fwrite(x = otu_css_norm, file = file.path(outdir, "otu_norm_CSS.csv"))

## relative abundances
otu_relative = transform_sample_counts(otu_tax_sample, function(x) x/sum(x) )
otu_rel_filtered = filter_taxa(otu_relative, function(x) mean(x) > 1e-3, TRUE)
nrow(otu_table(otu_rel_filtered))

writeLines(" - additional plots")
random_tree = rtree(ntaxa((otu_rel_filtered)), rooted=TRUE, tip.label=taxa_names(otu_rel_filtered))
plot(random_tree)

biom1 = merge_phyloseq(otu_rel_filtered, random_tree)
plot_tree(biom1, color=grouping_variable1, label.tips="taxa_names", ladderize="left", plot.margin=0.3)

png(filename = file.path(outdir, "figures", "genus_tree.png"), width = 1000, height = 800)
plot_tree(biom1, color="Genus", shape=grouping_variable1, size="abundance")
dev.off()

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
p <- plot_ordination(otu_tax_sample_norm, iMDS, color=grouping_variable1, shape=grouping_variable2)
ggsave(filename = file.path(outdir, "figures","mds_plot_bray_curtis.png"), plot = p, device = "png")

writeLines(" - write out distance matrix")
dd = dist2list(distances, tri = FALSE)
sample_order = gtools::mixedsort(as.character(unique(dd$row)), decreasing = TRUE)
dd$row = factor(dd$row, levels = sample_order)
dd$col = factor(dd$col, levels = sample_order)
dx <- dd |> pivot_wider(names_from = "col", values_from = "value")
fwrite(x = dx, file = file.path(outdir, "bray_curtis_distances.csv"))

## euclidean
writeLines(" - calculate Euclidean distances")
distances = distance(otu_tax_sample_norm, method="euclidean", type = "samples")
iMDS  <- ordinate(otu_tax_sample_norm, "MDS", distance=distances)
p <- plot_ordination(otu_tax_sample_norm, iMDS, color=grouping_variable1, shape=grouping_variable2)
ggsave(filename = file.path(outdir, "figures", "mds_plot_euclidean.png"), plot = p, device = "png")

writeLines(" - write out euclidean distance matrix")
dd = dist2list(distances, tri = FALSE)
dx <- dd |> pivot_wider(names_from = "col", values_from = "value")
fwrite(x = dx, file = file.path(outdir, "euclidean_distances.csv"))

### add tree ###
random_tree = rtree(ntaxa((otu_tax_sample_norm)), rooted=TRUE, tip.label=taxa_names(otu_tax_sample_norm))
otu_norm_tree = merge_phyloseq(otu_tax_sample_norm, random_tree)
# plot_tree(otu_norm_tree, color="treatment", label.tips="taxa_names", ladderize="left", plot.margin=0.3)
# plot_tree(otu_norm_tree, color="Phylum", shape=grouping_variable1, size="abundance")

## unifrac
writeLines(" - calculate Unifrac distances")
distances = distance(otu_norm_tree, method="unifrac", type = "samples")
iMDS  <- ordinate(otu_norm_tree, "MDS", distance=distances)
p <- plot_ordination(biom1, iMDS, color=grouping_variable1, shape=grouping_variable2)
# ggsave(filename = file.path(prj_folder, analysis_folder, "figures","mds_plot_unifrac.png"), plot = p, device = "png")

writeLines(" - write out Unifrac distance matrix")
dd = dist2list(distances, tri = FALSE)
dx <- dd |> pivot_wider(names_from = "col", values_from = "value")
fwrite(x = dx, file = file.path(outdir, "unifrac_distances.csv"))

## weighted unifrac
writeLines(" - calculate weighted Unifrac distances")
distances = distance(otu_norm_tree, method="wunifrac", type = "samples")
iMDS  <- ordinate(otu_norm_tree, "MDS", distance=distances)
p <- plot_ordination(otu_tax_sample_norm, iMDS, color=grouping_variable1, shape=grouping_variable2)
# ggsave(filename = file.path(prj_folder, analysis_folder, "figures","mds_plot_weighted_unifrac.png"), plot = p, device = "png")

writeLines(" - write out weighted Unifrac distance matrix")
dd = dist2list(distances, tri = FALSE)
dx <- dd |> pivot_wider(names_from = "col", values_from = "value")
fwrite(x = dx, file = file.path(outdir, "weighted_unifrac_distances.csv"))


