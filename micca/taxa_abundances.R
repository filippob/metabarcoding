
## SET UP
library("DT")
library("ape")
library("broom")
library("vegan")
library("ggpubr")
library("ggplot2")
library("ggridges")
library("tidytext")
library("phyloseq")
library("tidyverse")
library("data.table")
library("metagenomeSeq")

## PARAMETERS
## PARAMETERS
args = commandArgs(trailingOnly=TRUE)
if (length(args) >= 1) {
  
  #loading the parameters
  source(args[1])
  # source("Analysis/hrr/config.R")
  
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
    prjfolder = "Documents/cremonesi/suini_bontempo/pig_feces",
    analysis_folder = "Analysis/results",
    otu_norm_file = "otu_norm_CSS.csv",
    conf_file = "Config/rectum_mapping.csv",
    suffix = "feces_porous_zinc",
    nfactors = 2, ## n. of design variables (e.g. treatment and timpoint --> nfactors = 2)
    min_tot_n = 15,
    min_sample = 3,
    project = "ZnO poroso",
    sample_column = "sample",
    treatment_column = "treatment",
    grouping_variable2 = "timepoint",
    grouping_variable1 = "treatment",
    exp_levels = paste(c("CTR","T1","T2","T3"), collapse = ","),
    force_overwrite = FALSE
  ))
}

HOME <- Sys.getenv("HOME")
repo = file.path(HOME, config$repo)
prjfolder = file.path(HOME, config$prjfolder)
outdir = file.path(prjfolder,config$analysis_folder)
otu_norm_fname = file.path(prjfolder, config$analysis_folder, config$otu_norm_file)

fname = file.path(outdir, "taxa_abundance.config.r")
fwrite(x = config, file = fname)

## treatment levels as in the metadata file
print(config$exp_levels)
exp_levels = strsplit(config$exp_levels, split = ",")[[1]]

## read metadata
metadata = fread(file.path(prjfolder, config$conf_file))
metadata <- metadata |> rename(`sample-id` = !!config$sample_column, treatment = !!config$treatment_column)

if (config$project != "") metadata <- filter(metadata, project == !!config$project)

metadata <- mutate(metadata, `sample-id` = as.character(`sample-id`))
metadata = select(metadata, c(`sample-id`, treatment, timepoint))
metadata = metadata |> mutate(treatment = as.factor(treatment))

## read OTUs (normalised)
otus = fread(otu_norm_fname)
otus <- filter(otus, Phylum != "")

M <- as.data.frame(otus[,-c("tax_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")])

vec <- colnames(M) %in% paste("sample",metadata$`sample-id`,sep="-")
M <- M[,vec]

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

plot_title = paste("phylum", config$suffix)

phylums <- D |>
  group_by(Phylum) |>
  summarise(abundance = mean(avg)) |>
  arrange(desc(abundance))

D$Phylum <- factor(D$Phylum, levels = phylums$Phylum)

p <- ggplot(D, aes(x=factor(1), y=avg, fill=Phylum)) + geom_bar(width=1,stat="identity", alpha = 0.8)
p <- p + coord_polar(theta='y')
# p <- p + facet_wrap(~treatment)
p <- p + facet_grid(timepoint~treatment)
p <- p + xlab("relative abundances") + ylab("") + labs(title = plot_title)
p <- p + scale_fill_manual(values = my_palette)
p <- p + theme(text = element_text(size=12),
               axis.text.x = element_text(size=7),
               # axis.text.y = element_text(size=12),
               strip.text = element_text(size = 11),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               legend.text=element_text(size=8),
               legend.title=element_text(size=8),
               panel.spacing = unit(0.1, "lines")
)
print(p)

fname = paste("phylum_", config$suffix, ".png", sep="")
fname = file.path(outdir, "figures", fname)
ggsave(filename = fname, plot = p, device = "png", width = 7, height = 7)

fname = paste("phylum_relative_abundance_", config$suffix, ".csv", sep="")
fname = file.path(outdir, "tables", fname)
fwrite(x = arrange(D, timepoint, treatment, desc(avg)), file = fname, sep = ",")

to_save <- list("phylum_relabund"=D)

## F:B ratio
temp <- otus |>
  select(-c(tax_id,Kingdom, Class, Order, Family, Genus, Species)) |>
  gather(key = "sample", value = "counts", -c("Phylum")) |>
  filter(Phylum != "")

temp <- temp %>% inner_join(select(metadata, -`sample-id`), by = c("sample" = "sample"))

fb <- temp |>
  filter(Phylum %in% c("Firmicutes", "Bacteroidetes")) |>
  group_by(sample, Phylum) |>
  summarise(tot = sum(counts), treatment = unique(treatment), timepoint = unique(timepoint)) |>
  spread(key = Phylum, value = tot) |>
  mutate(FB = Firmicutes/Bacteroidetes, treatment = factor(treatment, levels = exp_levels))

to_save[["data_for_fb"]] = temp
to_save[["fb_res"]] = fb

p1 <- ggplot(filter(fb, FB != Inf), aes(x = FB, y = treatment, fill = treatment)) 
p1 <- p1 + geom_density_ridges(aes(point_color = treatment, point_fill = treatment, point_shape = treatment),
                               alpha = .2, point_alpha = 1, jittered_points = TRUE) 
p1 <- p1 + facet_wrap(~timepoint)
p1 <- p1 + scale_point_color_hue(l = 40)
p1 <- p1 + scale_discrete_manual(aesthetics = "point_shape", values = (seq(1,length(exp_levels))+20))
p1

fname = paste("fb_density_", config$suffix, ".png", sep="")
fname = file.path(outdir, fname)
ggsave(filename = fname, plot = p, device = "png")

temp <- filter(fb, FB != Inf) |>
  group_by(timepoint) |>
  do(tidy(lm(FB ~ treatment, data = .))) |>
  filter(term != "(Intercept)") |>
  mutate(term = gsub('treatment','',term)) |>
  rename(treatment = term) |>
  select(timepoint, treatment, p.value)

writeLines(" - FB ratio")
print(temp)

##############################
## taxonomic levels abundances
temp <- otus |>
  select(-c(tax_id,Kingdom, Phylum, Species)) |>
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


temp$level = factor(temp$level, levels = c("Class","Order","Family","Genus"))

temp <- temp |>
  mutate(taxon = reorder_within(taxon, avg_abund, level))

for (k in unique(temp$timepoint)) {
  
  print(paste("processing stratifying variable",k))
  
  p <- ggplot(filter(filter(temp, timepoint == k), avg_abund > 0.001), aes(x = taxon, y = avg_abund)) + geom_bar(aes(fill=taxon), stat = "identity")
  p <- p + coord_flip() + facet_grid(level~treatment, scales = "free", space = "free")
  p <- p + guides(fill="none") + scale_x_discrete(labels = function(x) str_replace(x, "_.*$", ""))
  p <- p + theme(axis.text.y = element_text(size = 4), axis.text.x = element_text(size = 5), axis.title.x = NULL)
  
  fname = paste("abundances_", config$suffix, "_", k, ".png", sep="")
  fname = file.path(outdir, "figures", fname)
  ggsave(filename = fname, plot = p, device = "png", width = 5, height = 9.5)
}


################
## linear models
################
otus <- filter(otus, Genus != "")
otus <- select(otus, -c(tax_id,Kingdom, Phylum, Class, Order, Family,Species))

mO <- otus |>
  gather(key = "sample", value = "counts", -c(Genus))

temp = select(metadata, c(sample, treatment, timepoint))
mO <- mO %>% inner_join(temp, by = c("sample" = "sample"))

## model y = mu + type + treatment + e
# mO$treatment <- factor(mO$treatment, levels = c("CTR", "T1", "T2", "T3"))
# mO$treatment <- factor(mO$treatment, levels = c("Controllo", "Lupini"))
mO$treatment <- factor(mO$treatment, levels = exp_levels)

dd <- mO |>
  group_by(Genus, treatment, timepoint) |>
  summarise(avg = mean(counts)) |>
  spread(key = "treatment", value = avg)

genus_stats <- mO %>% 
  group_by(Genus, treatment, timepoint) %>%
  summarise(avg = round(mean(counts),3), std = round(sd(counts),3))

diff <- genus_stats |> 
  select(-std) |> 
  spread(key = "treatment", value = avg)

cmb <- combn(exp_levels, 2)

res = data.frame("dummy" = rep("",nrow(diff)))
for(i in 1:ncol(cmb)) {
  
  pair <- cmb[,i]
  print(pair)
  
  x <- diff |>
    mutate(
      diff = .data[[pair[2]]] - .data[[pair[1]]]
    ) |> pull(diff)
  
  colname = paste(rev(pair), collapse="-")
  df = data.frame(a = x)
  names(df) <- colname
  
  res <- bind_cols(res,df)
}

res$dummy <- NULL

diff <- diff |> 
  bind_cols(res)

contrasts <- diff |>
  select(-all_of(exp_levels)) |>
  gather(key = "treatment", value = "difference", -c(Genus,timepoint))

to_save[["relabund_for_lm"]] = mO
to_save[["genus_stats"]] = genus_stats
to_save[["diffs"]] = contrasts

D <- mO %>%
  group_by(Genus,.data[[config$grouping_variable2]]) %>%
  do(tidy(lm(counts ~ .data[[config$grouping_variable1]], data=.))) %>%
  filter(term != "(Intercept)")

D$term <- gsub("treatment","",D$term)
D$term <- gsub("\\.data.*]]","",D$term)
# D$term <- gsub("timepoint","",D$term)

datatable(D, options = list(pageLength=100)) %>% 
  formatStyle('p.value', backgroundColor = styleInterval(0.05, c('yellow', 'white')))

fname = paste("significant_otus_treatment_within_timepoint_", config$suffix, ".csv", sep="")
fname = file.path(outdir, "tables", fname)
filter(D, p.value <= 0.05) %>% fwrite(file = fname, sep = ",", col.names = TRUE)

temp <- ungroup(D) |> filter(p.value <= 0.05)

tmp <- mO |> filter(Genus %in% temp$Genus) |>
  group_by(Genus, .data[[config$grouping_variable1]]) |>
  summarise(avg = mean(counts)) |>
  spread(key = .data[[config$grouping_variable1]], value = avg)

fname = paste("avg_counts_", config$suffix, ".csv", sep="")  
fname = file.path(outdir, "tables", fname)
fwrite(x = tmp, file = fname)

tmp <- mO |> filter(Genus %in% temp$Genus) |>
  group_by(Genus, .data[[config$grouping_variable2]], .data[[config$grouping_variable1]]) |>
  summarise(avg = mean(counts)) |>
  spread(key = .data[[config$grouping_variable1]], value = avg)

genus_stats <- genus_stats |>
  inner_join(temp, by = c("Genus" = "Genus","timepoint"="timepoint", "treatment"="term"))

fname = paste("significant_otus_abundance_", config$suffix, ".csv", sep="")
fname = file.path(outdir, "tables", fname)
fwrite(genus_stats, file = fname)

to_save[["lm_results"]] = D
to_save[["lm_significant_stats"]] = genus_stats

# ggplot(temp, aes(Genus)) + geom_bar()

# tmp <- temp |> group_by(Genus) |> summarise(N = n())
# 
# genus_stats <- genus_stats |> mutate(sign = ifelse(difference_vs_ctrl > 0, '+','-'))
# 
# g1 <- ggbarplot(genus_stats, x = "Genus", y = "difference_vs_ctrl", facet.by = "timepoint",
#           orientation = "horiz", fill="sign", color = "sign")
# 
# g1 <- ggpar(g1, font.tickslab = c(8))
# print(g1)
# 
# fname = paste("significant_otu_abundance_", config$suffix, ".png", sep="")
# fname = file.path(outdir, "figures", fname)
# ggsave(filename = fname, plot = g1, device = "png", width = 7.5, height = 7)
# 
# to_save[["significant_genus_stats"]] = genus_stats

########################################
## MANUALLY SET THE X VARIABLE IN AOV()
########################################
contrasts <- mO |>
  # nest(data = -c(timepoint,metric)) |>
  nest(data = -c(Genus, .data[[config$grouping_variable2]])) |>
  mutate(
    fit = map(data, ~ aov(counts ~ treatment, data = .x)),
    hsd = map(fit, TukeyHSD),
    tidied = map(hsd, tidy)
  ) |>
  unnest(tidied) |>
  select(-c(data,fit,hsd))

# vec <- contrasts$Genus %in% temp$Genus
# tmp <- contrasts[vec,]

to_save[["contrasts"]] = contrasts

fname = paste("taxa_contrasts_all_", config$suffix, ".csv", sep="")
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

fname = paste("taxa_contrasts_significant_", config$suffix, ".png", sep="")
fname = file.path(outdir, "figures", fname)
ggsave(filename = fname, plot = gg, device = "png", width = 5, height = 11)

## save results to R object
fname = paste("taxa_abundance_results_", config$suffix, ".RData", sep="")
fname = file.path(outdir, fname)
save(to_save, file = fname)

print("DONE!")
