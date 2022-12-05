
## SET UP
library("ape")
library("knitr")
library("vegan")
library("ggplot2")
library("phyloseq")
library("tidyverse")
library("data.table")
library("metagenomeSeq")

## PARAMETERS
HOME <- Sys.getenv("HOME")
prj_folder = file.path(HOME, "Documents/USEFUL")
analysis_folder = "Analysis/microbiome_chlorine/micca"
conf_file = "Config/mapping_file_chlorine.csv"
outdir = file.path(analysis_folder)

repo = file.path(prj_folder, "metabarcoding")
source(file.path(repo, "support_functions/dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(repo, "support_functions/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/

## loading data previously imported in phyloseq
load("Analysis/microbiome_chlorine/micca/otu_table/phyloseq.RData")

## Preprocessing: e.g. filtering
## making normalization folder
if(!file.exists(file.path(prj_folder, outdir))) dir.create(file.path(prj_folder, outdir), showWarnings = FALSE)

writeLines(" - CSS normalization")
otu_tax_sample_norm = phyloseq_transform_css(otu_tax_sample, norm = TRUE, log = FALSE)
# sample_data(otu_tax_sample_norm)
tipos = unique(sample_data(otu_tax_sample_norm)$type)

###############
## distances ##
###############
writeLines(" - beta diversity: distance matrices")
writeLines(" - available distance metrics")
## bray-curtis
for (k in tipos) {
  
  print(paste(" - calculate Bray-Curtis distances for ", k))
  temp = subset_samples(otu_tax_sample_norm, type == k)
  distances = distance(temp, method="bray", type = "samples")
  iMDS  <- ordinate(temp, "MDS", distance=distances)
  p <- plot_ordination(temp, iMDS, color="treatment")
  fname = paste("mds_plot_bray_curtis_", k , ".png")
  ggsave(filename = file.path(prj_folder, analysis_folder, "results", "figures", fname), plot = p, device = "png")
}

## permanova
metadata <- sample_data(otu_tax_sample_norm)
metadata$`sample-id` = row.names(metadata)
row.names(metadata) <- NULL
metadata <- as_tibble(metadata)

distances = distance(otu_tax_sample_norm, method="bray", type = "samples")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")

temp = dplyr::select(metadata, c(`sample-id`,type,treatment))
dx <- dx %>% inner_join(temp, by = c("row" = "sample-id"))

temp = select(dx,-row)
nvars = 2
matx= data.matrix(temp[,seq(1,ncol(temp)-nvars)])

print("Permanova on the whole dataset: comparison between types, treatmetns and type x treatment:")
source("~/Documents/cremonesi/rumine_anafi/parwise.adonis.r")
obj <- adonis2(matx ~ dx$type*dx$treatment, permutations = 1000)
print(kable(obj))

tipos = unique(dx$type)
for (tt in tipos) {
  
  print(tt)
  vec = (dx$type == tt)
  temp2 <- dx[vec,c(TRUE,vec,TRUE,TRUE)]
  temp1 <- matx[vec,vec]
  
  print(paste("Permanova on the", tt, "subset: comparison between treatments:"))
  obj <- adonis2(temp1 ~ temp2$treatment, permutations = 1000)
  print(kable(obj))
}

print("DONE!")
