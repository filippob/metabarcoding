
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
repo = file.path(HOME, "Documents/cremonesi/metabarcoding")
prj_folder = file.path(HOME, "Documents/cremonesi/suini_bontempo/pig_feces")
analysis_folder = "Analysis/results"
fname = "filtered_otu/otu_table_filtered.biom"
conf_file = "Config/rectum_mapping.csv"
outdir = file.path(analysis_folder)
nfactors = 2 ## n. of design variables (e.g. treatment and timpoint --> nfactors = 2)

source(file.path(repo, "support_functions/dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(repo, "support_functions/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/

## loading data previously imported in phyloseq
fname = file.path(prj_folder, analysis_folder, "phyloseq.RData")
load(fname)

## Preprocessing: e.g. filtering
## making normalization folder
if(!file.exists(file.path(prj_folder, outdir))) dir.create(file.path(prj_folder, outdir), showWarnings = FALSE)

writeLines(" - CSS normalization")
otu_tax_sample_norm = phyloseq_transform_css(otu_tax_sample, norm = TRUE, log = FALSE)
# sample_data(otu_tax_sample_norm)
tipos = unique(sample_data(otu_tax_sample_norm)$timepoint)

###############
## distances ##
###############
writeLines(" - beta diversity: distance matrices")
writeLines(" - available distance metrics")
## bray-curtis
for (k in tipos) {
  
  print(paste(" - calculate Bray-Curtis distances for ", k))
  temp = subset_samples(otu_tax_sample_norm, timepoint == k)
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

temp = dplyr::select(metadata, c(`sample-id`,timepoint,treatment))
dx <- dx %>% inner_join(temp, by = c("row" = "sample-id"))

temp = select(dx,-row)

print("Permanova on the whole dataset: comparison between types, treatmetns and type x treatment:")
source("~/Documents/cremonesi/rumine_anafi/parwise.adonis.r")

if (nfactors > 1) {
  
  nvars = nfactors
  matx= data.matrix(temp[,seq(1,ncol(temp)-nvars)])
  
  obj <- adonis2(matx ~ dx$timepoint*dx$treatment, permutations = 1000)
  print(kable(obj))
  
  tipos = unique(dx$timepoint)
  for (tt in tipos) {
    
    print(tt)
    vec = (dx$timepoint == tt)
    temp2 <- dx[vec,c(TRUE,vec,TRUE,TRUE)]
    temp1 <- matx[vec,vec]
    
    print(paste("Permanova on the", tt, "subset: comparison between treatments:"))
    obj <- adonis2(temp1 ~ temp2$treatment, permutations = 1000)
    print(kable(obj))
  }
} else {
  
  nvars = nfactors
  matx= data.matrix(temp[,seq(1,ncol(temp)-nvars)])
  
  obj <- adonis2(matx ~ dx$treatment, permutations = 1000)
  print(kable(obj))

  obj <- adonis2(matx ~ dx$treatment, permutations = 1000)
  print(kable(obj))
}



print("DONE!")
