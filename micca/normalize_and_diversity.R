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
HOME <- Sys.getenv("HOME")
prj_folder = file.path(HOME, "Documents/leguplus")
analysis_folder = "Analysis/results_lupini"
# conf_file = "Config/mapping_file.csv"
conf_file = "mapping_file.csv"
outdir = file.path(analysis_folder)

grouping_variable1 = "timepoint"
grouping_variable2 = "treatment"

repo = file.path(HOME, "Documents/cremonesi/metabarcoding")
source(file.path(repo, "support_functions/dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(repo, "support_functions/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/

## loading data previously imported in phyloseq
fname = file.path(prj_folder, analysis_folder, "phyloseq.RData")
load(fname)

## making results folder
# if(!file.exists(file.path(prj_folder, analysis_folder, "results"))) dir.create(file.path(prj_folder, analysis_folder, "results"), showWarnings = FALSE)

## making figures folder
if(!file.exists(file.path(prj_folder, analysis_folder, "figures"))) dir.create(file.path(prj_folder, analysis_folder, "figures"), showWarnings = FALSE)

## Alpha diversity
## alpha diversity is calculated on the original count data, not normalised 
## (see https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis)
## (see https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531)
writeLines(" - calculate alpha diversity indices")
alpha = estimate_richness(otu_tax_sample, split = TRUE)
alpha$"sample-id" = row.names(alpha)
alpha = relocate(alpha, `sample-id`)
fwrite(x = alpha, file = file.path(prj_folder, analysis_folder, "alpha.csv"))
p <- plot_richness(otu_tax_sample, x=grouping_variable2, color=grouping_variable1)
ggsave(filename = file.path(prj_folder, analysis_folder, "figures", "alpha_plot.png"), plot = p, device = "png", width = 11, height = 7)

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
fwrite(x = otu_css_norm, file = file.path(prj_folder, outdir, "otu_norm_CSS.csv"))

## relative abundances
otu_relative = transform_sample_counts(otu_tax_sample, function(x) x/sum(x) )
otu_rel_filtered = filter_taxa(otu_relative, function(x) mean(x) > 5e-3, TRUE)
nrow(otu_table(otu_rel_filtered))

writeLines(" - additional plots")
random_tree = rtree(ntaxa((otu_rel_filtered)), rooted=TRUE, tip.label=taxa_names(otu_rel_filtered))
plot(random_tree)

biom1 = merge_phyloseq(otu_rel_filtered, random_tree)
plot_tree(biom1, color=grouping_variable1, label.tips="taxa_names", ladderize="left", plot.margin=0.3)

png(filename = file.path(prj_folder, analysis_folder, "figures", "genus_tree.png"), width = 1000, height = 800)
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
ggsave(filename = file.path(prj_folder, analysis_folder, "figures","mds_plot_bray_curtis.png"), plot = p, device = "png")

writeLines(" - write out distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, analysis_folder, "bray_curtis_distances.csv"))

## euclidean
writeLines(" - calculate Euclidean distances")
distances = distance(otu_tax_sample_norm, method="euclidean", type = "samples")
iMDS  <- ordinate(otu_tax_sample_norm, "MDS", distance=distances)
p <- plot_ordination(otu_tax_sample_norm, iMDS, color=grouping_variable1, shape=grouping_variable2)
ggsave(filename = file.path(prj_folder, analysis_folder, "figures", "mds_plot_euclidean.png"), plot = p, device = "png")

writeLines(" - write out euclidean distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, analysis_folder, "euclidean_distances.csv"))

### add tree ###
random_tree = rtree(ntaxa((otu_tax_sample_norm)), rooted=TRUE, tip.label=taxa_names(otu_tax_sample_norm))
otu_norm_tree = merge_phyloseq(otu_tax_sample_norm, random_tree)
# plot_tree(otu_norm_tree, color="treatment", label.tips="taxa_names", ladderize="left", plot.margin=0.3)
plot_tree(otu_norm_tree, color="Phylum", shape=grouping_variable1, size="abundance")

## unifrac
writeLines(" - calculate Unifrac distances")
distances = distance(otu_norm_tree, method="unifrac", type = "samples")
iMDS  <- ordinate(otu_norm_tree, "MDS", distance=distances)
p <- plot_ordination(biom1, iMDS, color=grouping_variable1, shape=grouping_variable2)
ggsave(filename = file.path(prj_folder, analysis_folder, "figures","mds_plot_unifrac.png"), plot = p, device = "png")

writeLines(" - write out Unifrac distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, analysis_folder, "unifrac_distances.csv"))

## weighted unifrac
writeLines(" - calculate weighted Unifrac distances")
distances = distance(otu_norm_tree, method="wunifrac", type = "samples")
iMDS  <- ordinate(otu_norm_tree, "MDS", distance=distances)
p <- plot_ordination(biom1, iMDS, color=grouping_variable1, shape=grouping_variable2)
ggsave(filename = file.path(prj_folder, analysis_folder, "figures","mds_plot_weighted_unifrac.png"), plot = p, device = "png")

writeLines(" - write out weighted Unifrac distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(prj_folder, analysis_folder, "weighted_unifrac_distances.csv"))


