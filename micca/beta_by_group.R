
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
prj_folder = file.path(HOME, "Documents/cremonesi/suini_bontempo")
analysis_folder = "Analysis/micca/results_zinc_caecum"
fname = "filtered_otu/otu_table_filtered.biom"
conf_file = "Config/zinco_poroso_caecum.csv"
outdir = file.path(prj_folder,analysis_folder)
suffix = "zinc_caecum"
nfactors = 1 ## n. of design variables (e.g. treatment and timpoint --> nfactors = 2)

source(file.path(repo, "support_functions/dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(repo, "support_functions/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/

## loading data previously imported in phyloseq
fname = file.path(outdir, "phyloseq.RData")
load(fname)

## Preprocessing: e.g. filtering
## making normalization folder
if(!file.exists(file.path(outdir))) dir.create(file.path(outdir), showWarnings = FALSE)

writeLines(" - CSS normalization")
otu_tax_sample_norm = phyloseq_transform_css(otu_tax_sample, norm = TRUE, log = FALSE)
# sample_data(otu_tax_sample_norm)

## subset samples based on metadata
# otu_tax_sample_norm = subset_samples(otu_tax_sample_norm, project == "ZnO poroso")

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
  ggsave(filename = file.path(prj_folder, analysis_folder, "figures", fname), plot = p, device = "png")
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

fpath = file.path(outdir, "tables")
dir.create(fpath, showWarnings = FALSE)
fname = file.path(fpath, paste("permanova_",suffix,".csv",sep=""))
fwrite(x = list("PERMANOVA"), file = fname)

write(x = "PERMANOVA", file = fname)

if (nfactors > 1) {
  
  nvars = nfactors
  matx= data.matrix(temp[,seq(1,ncol(temp)-nvars)])
  
  obj <- adonis2(matx ~ dx$timepoint+dx$treatment, permutations = 1000)
  print(kable(obj))
  # fwrite(x = list(kable(obj)), file = fname, append = TRUE)
  write(x = "ALL DATA", file = fname, append = TRUE)
  write(x = kable(obj), file = fname, append = TRUE)
  
  tipos = unique(dx$timepoint)
  for (tt in tipos) {
    
    print(tt)
    vec = (dx$timepoint == tt)
    temp2 <- dx[vec,c(TRUE,vec,TRUE,TRUE)]
    temp1 <- matx[vec,vec]
    
    print(paste("Permanova on the", tt, "subset: comparison between treatments:"))
    obj <- adonis2(temp1 ~ temp2$treatment, permutations = 1000)
    write(x = tt, file = fname, append = TRUE)
    write(x = kable(obj), file = fname, append = TRUE)
    print(kable(obj))
  }
} else {
  
  nvars = nfactors
  matx= data.matrix(temp[,seq(1,ncol(temp)-nvars)])
  
  obj <- adonis2(matx ~ dx$treatment, permutations = 1000)
  print(kable(obj))
  write(x = kable(obj), file = fname, append = TRUE)
}

library("ggpubr")
library("ggfortify")

for (k in tipos) {
  
  print(paste("processing stratifying variable",k))
  temp <- subset_samples(otu_tax_sample_norm, timepoint == k)
  mtd <- sample_data(temp)
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
                 repel = TRUE)
  
  fname = paste("mds_plot_bray-curtis_ellipse_", suffix, "_", k , ".png")
  ggsave(filename = file.path(prj_folder, analysis_folder, "figures", fname), plot = g, device = "png")
}




g

print("DONE!")
