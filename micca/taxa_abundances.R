
## SET UP
library("ape")
library("broom")
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

## read metadata
metadata = fread(file.path(prj_folder, conf_file))

## read OTUs (normalised)
fname = file.path(prj_folder, analysis_folder, "results/otu_norm_CSS.csv")
otus = fread(fname)
otus <- filter(otus, Phylum != "")

M <- otus[,-c("tax_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
## !! remember: due to how matrices are stored internally in R, M/colSums(M) won't give the expected results --> use sweep instead !!
M <- sweep(M, MARGIN=2, FUN="/", STATS=colSums(M))
M <- cbind.data.frame("Phylum"=otus$Phylum, M)

mO <- M %>% gather(key = "sample", value = "counts", -c("Phylum"))

metadata$sample <- paste("sample", metadata$`sample-id`, sep = "-")
temp = select(metadata, c(sample, type, treatment))
mO <- mO %>% inner_join(temp, by = "sample")

D <- mO %>% 
  group_by(type,treatment,Phylum) %>% 
  summarise(avg = mean(counts))


library("ggpubr")
my_palette = get_palette(c("#00AFBB", "#E7B800", "#FC4E07","darkgrey"), length(unique(D$Phylum)))

p <- ggplot(D, aes(x=factor(1), y=avg, fill=Phylum)) + geom_bar(width=1,stat="identity")
p <- p + coord_polar(theta='y')
p <- p + facet_grid(type~treatment)
p <- p + xlab("relative abundances") + ylab("") + labs(title = "cow's milk")
p <- p + scale_fill_manual(values = my_palette)
p <- p + theme(text = element_text(size=12),
               axis.text.x = element_text(size=6),
               # axis.text.y = element_text(size=12),
               strip.text = element_text(size = 12),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               legend.text=element_text(size=7),
               legend.title=element_text(size=8),
               panel.spacing = unit(0.1, "lines")
)

print(p)

## linear model
fname = file.path(prj_folder, analysis_folder, "results/otu_norm_CSS.csv")
otus = fread(fname)
otus <- filter(otus, Genus != "")
otus <- select(otus, -c(tax_id,Kingdom, Phylum, Class, Order, Family, Species))

mO <- otus |>
  gather(key = "sample", value = "counts", -Genus)

temp = select(metadata, c(sample, type, treatment))
mO <- mO %>% inner_join(temp, by = c("sample" = "sample"))

## model y = mu + type + treatment + e
mO$treatment <- factor(mO$treatment, levels = c("NO CLORO", "CLORO"))

temp = filter(mO, Genus == "Methylobacterium")
fit <- lm(counts ~ type + treatment, data = temp)
glance(fit)
tidy(fit)
tidy(anova(fit))

D <- mO %>%
  group_by(Genus) %>%
  do(tidy(lm(counts ~  type + treatment, data=.))) %>%
  # do(tidy(anova(lm(counts ~  type + treatment, data=.)))) %>%
  filter(term == "treatmentCLORO")

library("DT")
datatable(D, options = list(pageLength=100)) %>% 
  formatStyle('p.value', backgroundColor = styleInterval(0.05, c('yellow', 'white')))

fname = file.path(prj_folder, analysis_folder, "significant_otus_treatment_plus_type.csv")
filter(D, p.value <= 0.05) %>%fwrite(file = fname, sep = ",", col.names = TRUE)

## model y = mu + treatment + e (within type)
D <- mO %>%
  group_by(type,Genus) %>%
  do(tidy(lm(counts ~ treatment, data=.))) %>%
  filter(term == "treatmentCLORO")

siero = filter(D, type == "SIEROINNESTO")
fwrite(x = siero, file = "Analysis/microbiome_chlorine/micca/results/siero_innesto.csv")

fname = file.path(prj_folder, analysis_folder, "significant_otus_treatment_within_type.csv")
filter(D, p.value <= 0.05) %>%fwrite(file = fname, sep = ",", col.names = TRUE)
