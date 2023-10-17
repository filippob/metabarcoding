
## SET UP
library("ape")
library("broom")
library("vegan")
library("ggpubr")
library("ggplot2")
library("phyloseq")
library("tidyverse")
library("data.table")
library("metagenomeSeq")

## PARAMETERS
HOME <- Sys.getenv("HOME")
repo = file.path(HOME, "Documents/cremonesi/metabarcoding")
prj_folder = file.path(HOME, "Documents/leguplus")
analysis_folder = "Analysis/results_lupini"
otu_norm_fname = file.path(prj_folder, analysis_folder, "otu_norm_CSS.csv")
conf_file = "mapping_file.csv"
outdir = file.path(prj_folder,analysis_folder)
suffix = "lupino"

## treatment levels as in the metadata file
level1 = "Controllo"
level2 = "Lupini"

## read metadata
metadata = fread(file.path(prj_folder, conf_file))
# metadata <- metadata |> select(`progressivo MISEQ`, `GRUPPO sperimentale`, TIMEPOINT, Progetto, tipologia)
# metadata <- rename(metadata, `sample-id` = `progressivo MISEQ`, treatment = `GRUPPO sperimentale`, timepoint = TIMEPOINT)
metadata <- metadata |> mutate(`sample-id` = as.character(`sample-id`))
# metadata <- metadata |> rename(`sample-id` = id)
# metadata <- mutate(metadata, id = as.character(id))
metadata = select(metadata, c(`sample-id`, timepoint, treatment))
# metadata = select(metadata, c(`sample-id`, treatment, timepoint)) |> filter(timepoint == "D28")

## read OTUs (normalised)
otus = fread(otu_norm_fname)
otus <- filter(otus, Phylum != "")

# M <- otus[,-c("tax_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
M <- otus[,-c("tax_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")]
## !! remember: due to how matrices are stored internally in R, M/colSums(M) won't give the expected results --> use sweep instead !!
M <- sweep(M, MARGIN=2, FUN="/", STATS=colSums(M))
M <- cbind.data.frame("Phylum"=otus$Phylum, M)

metadata$sample <- paste("sample", metadata$`sample-id`, sep = "-")
mO <- M %>% gather(key = "sample", value = "counts", -c("Phylum"))

temp = select(metadata, c(sample, treatment, timepoint))
mO <- mO %>% inner_join(temp, by = "sample")

D <- mO %>% 
  group_by(treatment,timepoint,Phylum) %>% 
  dplyr::summarise(avg = sum(counts))

D <- mO |> group_by(treatment, timepoint) |>
  mutate(tot = sum(counts)) |>
  group_by(treatment,timepoint, Phylum) |>
  summarise(avg = sum(counts)/unique(tot))

  
my_palette = get_palette(c("#00AFBB", "pink", "#FC4E07", "green", "purple", "#E7B800", "darkgrey", "salmon", "blue"), length(unique(D$Phylum)))

plot_title = paste("phylum", suffix)

p <- ggplot(D, aes(x=factor(1), y=avg, fill=Phylum)) + geom_bar(width=1,stat="identity", alpha = 0.8)
p <- p + coord_polar(theta='y')
# p <- p + facet_wrap(~treatment)
p <- p + facet_grid(timepoint~treatment)
p <- p + xlab("relative abundances") + ylab("") + labs(title = plot_title)
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

fname = paste("phylum_", suffix, ".png", sep="")
fname = file.path(outdir, "figures", fname)
ggsave(filename = fname, plot = p, device = "png", width = 7, height = 7)

fname = paste("phylum_relative_abundance_", suffix, ".csv", sep="")
fname = file.path(outdir, "tables", fname)
fwrite(x = D, file = fname, sep = ",")

## taxonomic levels abundances
otus = fread(otu_norm_fname)
temp <- otus |>
  select(-c(tax_id,Kingdom, Phylum)) |>
  gather(key = "sample", value = "counts", -c("Class","Order","Family","Genus")) |>
  gather(key = "level", value = "taxon", -c("sample","counts")) |>
  filter(taxon != "")

temp <- temp %>% inner_join(select(metadata, -`sample-id`), by = c("sample" = "sample"))

temp <- temp |> 
  group_by(sample, timepoint, treatment, level, taxon) |> 
  summarise(tot_sample = sum(counts)) |> 
  mutate(tot = sum(tot_sample), abundance = tot_sample/tot) |>
  group_by(timepoint, treatment, level, taxon) |>
  summarise(avg_abund = mean(abundance))

# otus |>
#   select(-c(tax_id,Kingdom, Phylum, Species, Order, Family, Genus)) |>
#   gather(key = "sample", value = "counts", -c("Class")) |>
#   group_by(sample, Class) |> 
#   summarise(tot_sample = sum(counts)) |> 
#   filter(Class != "") |>
#   mutate(tot = sum(tot_sample), abundance = tot_sample/tot) |>
#   group_by(Class) |>
#   summarise(avg_abund = mean(abundance))

temp$level = factor(temp$level, levels = c("Class","Order","Family","Genus"))

library("tidytext")
temp <- temp |>
  mutate(taxon = reorder_within(taxon, avg_abund, level))

for (k in unique(temp$timepoint)) {
  
  print(paste("processing stratifying variable",k))
  
  p <- ggplot(filter(filter(temp, timepoint == k), avg_abund > 0.001), aes(x = taxon, y = avg_abund)) + geom_bar(aes(fill=taxon), stat = "identity")
  p <- p + coord_flip() + facet_grid(level~treatment, scales = "free", space = "free")
  p <- p + guides(fill="none") + scale_x_discrete(labels = function(x) str_replace(x, "_.*$", ""))
  p <- p + theme(axis.text.y = element_text(size = 4), axis.text.x = element_text(size = 5), axis.title.x = NULL)
  
  fname = paste("abundances_", suffix, "_", k, ".png", sep="")
  fname = file.path(outdir, "figures", fname)
  ggsave(filename = fname, plot = p, device = "png", width = 5, height = 9.5)
}


################
## linear models
################
otus <- filter(otus, Genus != "")
# otus <- select(otus, -c(tax_id,Kingdom, Phylum, Class, Order, Family, Species))
otus <- select(otus, -c(tax_id,Kingdom, Phylum, Class, Order, Family))

mO <- otus |>
  gather(key = "sample", value = "counts", -c(Genus))

## We consider that the species level is not usable (virtually all are uncultered or unidentified)
# table(mO$Species)
# mO$Species <- NULL

temp = select(metadata, c(sample, treatment, timepoint))
mO <- mO %>% inner_join(temp, by = c("sample" = "sample"))

## model y = mu + type + treatment + e
# mO$treatment <- factor(mO$treatment, levels = c("CTR", "T1", "T2", "T3"))
# mO$treatment <- factor(mO$treatment, levels = c("Controllo", "Lupini"))
mO$treatment <- factor(mO$treatment, levels = c(level1, level2))

dd <- mO |>
  group_by(Genus, treatment, timepoint) |>
  summarise(avg = mean(counts)) |>
  spread(key = "treatment", value = avg)

temp = filter(mO, Genus == mO$Genus[1])
# fit <- lm(counts ~ timepoint + treatment, data = temp)
# glance(fit)
# tidy(fit)
# tidy(anova(fit))

genus_stats <- mO %>% 
  group_by(Genus, treatment, timepoint) %>%
  summarise(avg = round(mean(counts),3), std = round(sd(counts),3))

# temp <- genus_stats |> select(-std) |> spread(key = "treatment", value = avg) |>
#   mutate(diff = INSETTI-SOIA) |>
#   select(-c(SOIA,INSETTI))

# temp <- genus_stats |> select(-std) |> spread(key = "treatment", value = avg) |>
#   mutate(diffT1 = T1-CTR, diffT2 = T2-CTR, diffT3 = T3-CTR) |>
#   select(-c(CTR,T1,T2,T3)) |>
#   rename(T1 = diffT1, T2 = diffT2, T3 = diffT3) |>
#   gather(key = "treatment", value = "difference_vs_ctrl", -c(Genus,timepoint))

temp <- genus_stats |> select(-std) |> spread(key = "treatment", value = avg) |>
  mutate(diff = .data[[level2]]-.data[[level1]]) |>
  select(-all_of(c(level1, level2))) |>
  rename(!!level2 := diff) |>
  gather(key = "treatment", value = "difference_vs_ctrl", -c(Genus,timepoint))

genus_stats <- temp |>
  inner_join(genus_stats, by = c("Genus" = "Genus", "timepoint" = "timepoint", "treatment" = "treatment")) |>
  # inner_join(genus_stats, by = c("Genus" = "Genus", "treatment" = "treatment")) |>
  select(-avg)

D <- mO %>%
  group_by(Genus,timepoint) %>%
  do(tidy(lm(counts ~ treatment, data=.))) %>%
  # do(tidy(lm(counts ~ timepoint + treatment, data=.))) %>%
  # do(tidy(anova(lm(counts ~  type + treatment, data=.)))) %>%
  filter(term != "(Intercept)")

D$term <- gsub("treatment","",D$term)
# D$term <- gsub("timepoint","",D$term)

library("DT")
datatable(D, options = list(pageLength=100)) %>% 
  formatStyle('p.value', backgroundColor = styleInterval(0.05, c('yellow', 'white')))

fname = paste("significant_otus_treatment_within_timepoint_", suffix, ".csv", sep="")
fname = file.path(outdir, "tables", fname)
filter(D, p.value <= 0.05) %>% fwrite(file = fname, sep = ",", col.names = TRUE)

temp <- ungroup(D) |> filter(p.value <= 0.05)

tmp <- mO |> filter(Genus %in% temp$Genus) |>
  group_by(Genus, treatment) |>
  summarise(avg = mean(counts)) |>
  spread(key = treatment, value = avg)

fname = paste("avg_counts_", suffix, ".csv", sep="")  
fname = file.path(outdir, "tables", fname)
fwrite(x = tmp, file = fname)

tmp <- mO |> filter(Genus %in% temp$Genus) |>
  group_by(Genus, timepoint, treatment) |>
  summarise(avg = mean(counts)) |>
  spread(key = treatment, value = avg)

genus_stats <- genus_stats |>
  inner_join(tmp, by = c("Genus","timepoint"))

fname = paste("significant_otus_abundance_", suffix, ".csv", sep="")
fname = file.path(outdir, "tables", fname)
fwrite(genus_stats, file = fname)

ggplot(temp, aes(Genus)) + geom_bar()

library("ggpubr")
tmp <- temp |> group_by(Genus) |> summarise(N = n())
ggbarplot(tmp, "Genus", "N", sort.val = "desc")

?ggbarplot

genus_stats <- genus_stats |> mutate(sign = ifelse(difference_vs_ctrl > 0, '+','-'))

g1 <- ggbarplot(genus_stats, x = "Genus", y = "difference_vs_ctrl", facet.by = "timepoint",
          orientation = "horiz", fill="sign", color = "sign")

g1 <- ggpar(g1, font.tickslab = c(8))
print(g1)

fname = paste("significant_otu_abundance_", suffix, ".png", sep="")
fname = file.path(outdir, "figures", fname)
ggsave(filename = fname, plot = g1, device = "png", width = 7.5, height = 7)

## model y = mu + treatment + e (within type)
# D <- mO %>%
#   group_by(Genus) %>%
#   do(tidy(lm(counts ~ treatment, data=.))) %>%
#   filter(term == "(Intercept)")
# 
# 
# fname = file.path(prj_folder, analysis_folder, "significant_otus_treatment_within_type.csv")
# filter(D, p.value <= 0.05) %>%fwrite(file = fname, sep = ",", col.names = TRUE)

########################################
## MANUALLY SET THE X VARIABLE IN AOV()
########################################
contrasts <- mO |>
  # nest(data = -c(timepoint,metric)) |>
  nest(data = -c(Genus, timepoint)) |>
  mutate(
    fit = map(data, ~ aov(counts ~ treatment, data = .x)),
    hsd = map(fit, TukeyHSD),
    tidied = map(hsd, tidy)
  ) |>
  unnest(tidied) |>
  select(-c(data,fit,hsd))

# vec <- contrasts$Genus %in% temp$Genus
# tmp <- contrasts[vec,]

fname = paste("taxa_contrasts_all_", suffix, ".csv", sep="")
fname = file.path(outdir, "tables", fname)
fwrite(contrasts, file = fname)

tmp <- contrasts |>
  filter(!is.na(adj.p.value)) |>
  select(Genus, timepoint, contrast, adj.p.value) |>
  mutate(pvalue = ifelse(adj.p.value <= 0.05, "<= 0.05", "> 0.05")) |>
  select(-adj.p.value)

colors = c("red","white")

gg <- ggplot(tmp, aes(x = contrast, y = Genus)) + geom_tile(aes(fill=pvalue), color = "white") + 
  facet_wrap(~timepoint) + 
  # scale_fill_distiller(palette = "YlGnBu") + 
  scale_fill_manual(values=colors) + 
  theme(axis.text.x = element_text(angle=90, size = 6),
        axis.text.y = element_text(size = 4),
        legend.key.size = unit(0.1, 'cm'),
        strip.text = element_text(size=7),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6),
        legend.key = element_rect(color="black"),
        panel.background = element_rect(fill = 'white', color = 'white'),
        panel.grid.major = element_line(color = 'white'),
        panel.grid.minor = element_line(color = 'white'))

fname = paste("taxa_contrasts_significant_", suffix, ".png", sep="")
fname = file.path(outdir, "figures", fname)
ggsave(filename = fname, plot = gg, device = "png", width = 5, height = 11)
