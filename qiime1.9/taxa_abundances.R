
## SET UP
library("ape")
library("vegan")
library("ggplot2")
library("phyloseq")
library("tidyverse")
library("data.table")
library("metagenomeSeq")

## PARAMETERS
HOME <- Sys.getenv("HOME")
prj_folder = file.path(HOME, "Documents/USEFUL")
analysis_folder = "Analysis/microbiome_predipping/qiime1.9/results"
conf_file = "Config/mapping_file.csv"
outdir = file.path(analysis_folder)

repo = file.path(prj_folder, "metabarcoding")

## read metadata
metadata = fread(file.path(prj_folder, conf_file))
metadata <- rename(metadata, sample = `#sampleID`)
metadata <- metadata |>
  mutate(n = nchar(sample), numid = ifelse(n==2, paste("sample",substr(sample,1,1),sep="-"), paste("sample",substr(sample,3,4),sep="-")))

## read OTUs (normalised)
fname = file.path(prj_folder, analysis_folder, "otu_norm_CSS.csv")
otus = fread(fname)
otus <- filter(otus, Genus != "")
otus <- select(otus, -c(tax_id,Kingdom, Phylum, Class, Order, Family, Species))

mO <- otus |>
  gather(key = "sample", value = "counts", -Genus)

temp = select(metadata, c(numid, azienda, treatment))
mO <- mO %>% inner_join(temp, by = c("sample" = "numid"))

mO$treatment <- factor(mO$treatment, levels = c("no-cloro", "cloro"))

D <- mO %>%
  group_by(Genus) %>%
  do(tidy(lm(counts ~  azienda + treatment, data=.))) %>%
  filter(term == "treatmentcloro")

library("DT")
datatable(D, options = list(pageLength=100)) %>% 
  formatStyle('p.value', backgroundColor = styleInterval(0.05, c('yellow', 'white')))

fname = file.path(prj_folder, analysis_folder, "significant_otus_treatment_plus_herd.csv")
filter(D, p.value <= 0.05) %>%fwrite(file = fname, sep = ",", col.names = TRUE)


D <- mO %>%
  group_by(azienda,Genus) %>%
  do(tidy(lm(counts ~ treatment, data=.))) %>%
  filter(term == "treatmentcloro")

fname = file.path(prj_folder, analysis_folder, "significant_otus_treatment_within_herd(effect).csv")
filter(D, p.value <= 0.05) %>%fwrite(file = fname, sep = ",", col.names = TRUE)
