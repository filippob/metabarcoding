
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
    prjfolder = "Documents/cremonesi/sicurprolat",
    analysis_folder = "Analysis",
    otu_norm_file = "otu_norm_CSS.csv",
    conf_file = "Config/mapping_file.csv",
    suffix = "bruna",
    nfactors = 2, ## n. of design variables (e.g. treatment and timpoint --> nfactors = 2)
    min_tot_n = 15,
    min_sample = 3,
    project = "",
    sample_column = "sample-id",
    sample_prefix = "sample-",
    treatment_column = "season",
    grouping_variable2 = "place",
    grouping_variable1 = "treatment",
    exp_levels = paste(c("Freddo", "Caldo"), collapse = ","), ## !! THE FIRST LEVEL IS THE BENCHMARK !! Not treated,Treated
    sig_threshold = 0.05,
    force_overwrite = FALSE
  ))
}

HOME <- Sys.getenv("HOME")
repo = file.path(HOME, config$repo)
prjfolder = file.path(HOME, config$prjfolder)
outdir = file.path(prjfolder,config$analysis_folder)
otu_norm_fname = file.path(prjfolder, config$analysis_folder, config$otu_norm_file)

fname = file.path(outdir, "taxa_abundance.config.RData")
save(config, file = fname)

## treatment levels as in the metadata file
print(config$exp_levels)
exp_levels = strsplit(config$exp_levels, split = ",")[[1]]

## read metadata
metadata = fread(file.path(prjfolder, config$conf_file))
# if(config$treatment_column != "treatment" & "treatment" %in% names(metadata)) metadata <- rename(metadata, secondary_treatment = treatment)
if(config$treatment_column != "treatment" & "treatment" %in% names(metadata)) metadata$treatment <- NULL
metadata <- metadata |> rename(`sample-id` = !!config$sample_column, treatment = !!config$treatment_column)

if (config$project != "") metadata <- filter(metadata, project == !!config$project)

metadata <- mutate(metadata, `sample-id` = as.character(`sample-id`))
if(config$grouping_variable2 %in% names(metadata)) {
  metadata = metadata |> select(c(`sample-id`, treatment, all_of(config$grouping_variable2))) |> mutate(treatment = as.factor(treatment))
} else metadata = metadata |> select(c(`sample-id`, treatment, !!config$grouping_variable2))

## convert to factors
metadata <- metadata |> mutate_at(vars(-("sample-id")),as.factor)

## read OTUs (normalised)
otus = fread(otu_norm_fname)
otus <- filter(otus, Phylum != "")

M <- as.data.frame(otus[,-c("tax_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")])

vec <- colnames(M) %in% paste(config$sample_prefix,metadata$`sample-id`,sep="")
M <- M[,vec]

## !! remember: due to how matrices are stored internally in R, M/colSums(M) won't give the expected results --> use sweep instead !!
M <- sweep(M, MARGIN=2, FUN="/", STATS=colSums(M))
M <- cbind.data.frame("Phylum"=otus$Phylum, M)

metadata$sample <- paste(config$sample_prefix, metadata$`sample-id`, sep = "")
mO <- M %>% gather(key = "sample", value = "counts", -c("Phylum"))

if(config$grouping_variable2 %in% names(metadata)) {
  temp = select(metadata, c(sample, treatment, all_of(config$grouping_variable2)))
} else temp = select(metadata, c(sample, treatment, !!config$grouping_variable2))

mO <- mO %>% inner_join(temp, by = "sample")

if(config$grouping_variable2 %in% names(metadata)) {
  
  D <- mO %>% 
    group_by(treatment, .data[[config$grouping_variable2]], Phylum) %>% 
    dplyr::summarise(avg = sum(counts))
  
  D <- mO |> group_by(treatment, .data[[config$grouping_variable2]]) |>
    mutate(tot = sum(counts)) |>
    group_by(treatment, .data[[config$grouping_variable2]], Phylum) |>
    summarise(avg = sum(counts)/unique(tot))
} else {
  
  D <- mO %>% 
    dplyr::group_by(treatment, .data[[config$grouping_variable2]], Phylum) %>% 
    dplyr::summarise(avg = sum(counts))
  
  D <- mO |> group_by(treatment, .data[[config$grouping_variable2]]) |>
    mutate(tot = sum(counts)) |>
    group_by(treatment, .data[[config$grouping_variable2]], Phylum) |>
    summarise(avg = sum(counts)/unique(tot))
}

  
my_palette = get_palette(c("#00AFBB", "pink", "#FC4E07", "green", "purple", "#E7B800", "darkgrey", "salmon", "blue"), length(unique(D$Phylum)))

plot_title = paste("phylum", config$suffix)

phylums <- D |>
  group_by(Phylum) |>
  summarise(abundance = mean(avg)) |>
  arrange(desc(abundance))

D$Phylum <- factor(D$Phylum, levels = phylums$Phylum)

var1 = config$grouping_variable1
var2 = config$grouping_variable2

p <- ggplot(D, aes(x=factor(1), y=avg, fill=Phylum)) + geom_bar(width=1,stat="identity", alpha = 0.8)
p <- p + coord_polar(theta='y')
if (config$grouping_variable2 != "") {
  p <- p + facet_grid(rows = vars(!!sym(var1)), cols = vars(!!sym(var2)))
} else p <- p + facet_wrap(~treatment)
p <- p + xlab("relative abundances") + ylab(NULL) + labs(title = plot_title)
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
fwrite(x = arrange(D, treatment, desc(avg)), file = fname, sep = ",")

to_save <- list("phylum_relabund"=D)


##############################
## taxonomic levels abundances
temp <- otus |>
  select(-c(tax_id,Kingdom, Phylum, Species)) |>
  gather(key = "sample", value = "counts", -c("Class","Order","Family","Genus")) |>
  gather(key = "level", value = "taxon", -c("sample","counts")) |>
  filter(taxon != "")

temp <- temp %>% inner_join(select(metadata, -`sample-id`), by = c("sample" = "sample"))

if (config$grouping_variable2 != "") {
  temp <- temp |> 
    group_by(sample, .data[[config$grouping_variable2]], treatment, level, taxon) |> 
    summarise(tot_sample = sum(counts)) |> 
    mutate(tot = sum(tot_sample), abundance = tot_sample/tot) |>
    group_by(.data[[config$grouping_variable2]], treatment, level, taxon) |>
    summarise(avg_abund = mean(abundance), std = sd(abundance))
  
  
  temp$level = factor(temp$level, levels = c("Class","Order","Family","Genus"))
  
  temp <- temp |>
    mutate(taxon = reorder_within(taxon, avg_abund, level))
  
  steps = pull(temp, .data[[config$grouping_variable2]])
  for (k in unique(steps)) {
    
    print(paste("processing stratifying variable",k))
    
    p <- ggplot(filter(filter(temp, .data[[config$grouping_variable2]] == k), avg_abund > 0.001), aes(x = taxon, y = avg_abund)) + geom_bar(aes(fill=taxon), stat = "identity")
    p <- p + coord_flip() + facet_grid(level~treatment, scales = "free", space = "free")
    p <- p + guides(fill="none") + scale_x_discrete(labels = function(x) str_replace(x, "_.*$", ""))
    p <- p + theme(axis.text.y = element_text(size = 4), axis.text.x = element_text(size = 5), axis.title.x = NULL)
    
    fname = paste("abundances_", config$suffix, "_", k, ".png", sep="")
    fname = file.path(outdir, "figures", fname)
    ggsave(filename = fname, plot = p, device = "png", width = 5, height = 9.5)
  }
}

################
## linear models
################
otus <- filter(otus, Genus != "")
otus <- select(otus, -c(tax_id,Kingdom, Phylum, Class, Order, Family,Species))

mO <- otus |>
  gather(key = "sample", value = "counts", -c(Genus))

if (config$grouping_variable2 != "") { 
  temp = select(metadata, c(sample, treatment, all_of(config$grouping_variable2)))
} else temp = select(metadata, c(sample, treatment))
mO <- mO %>% inner_join(temp, by = c("sample" = "sample"))
mO$treatment = factor(mO$treatment, levels = exp_levels)

if (config$grouping_variable2 != "") { 

  temp <- mO |>
    group_by(Genus, treatment, .data[[config$grouping_variable2]]) |>
    summarise(avg = mean(counts))
  
  dd <- pivot_wider(temp, names_from = treatment, values_from = avg)
  
  genus_stats <- mO %>% 
    group_by(Genus, treatment, .data[[config$grouping_variable2]]) %>%
    summarise(avg = round(mean(counts),3), std = round(sd(counts),3), N=n())
} else {
  
  dd <- mO |>
    group_by(Genus, treatment) |>
    summarise(avg = mean(counts)) |>
    spread(key = "treatment", value = avg)
  
  genus_stats <- mO %>% 
    group_by(Genus, treatment) %>%
    summarise(avg = round(mean(counts),3), std = round(sd(counts),3), N = n())
}

diff <- genus_stats |> 
  select(-c(std,N)) |> 
  pivot_wider(names_from = treatment, values_from = avg)

std <- genus_stats |> 
  select(-c(avg,N)) |> 
  pivot_wider(names_from = treatment, values_from = std)

sample_size <- genus_stats |> 
  select(-c(avg,std)) |> 
  pivot_wider(names_from = treatment, values_from = N)

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

if (config$grouping_variable2 != "") { 
  contrasts <- diff |>
    select(-all_of(exp_levels)) |>
    gather(key = "treatment", value = "difference", -c(Genus, .data[[config$grouping_variable2]]))
} else {
  contrasts <- diff |>
    select(-all_of(exp_levels)) |>
    gather(key = "treatment", value = "difference", -c(Genus))
}


contrasts <- contrasts |> inner_join(std, by = c("Genus", config$grouping_variable2)) |>
  rename_with(~str_c("std_", .), .cols = all_of(exp_levels))
  
 contrasts <- contrasts |> inner_join(sample_size, by = c("Genus", config$grouping_variable2)) |>
   rename_with(~str_c("n_", .), .cols = all_of(exp_levels))
  # rename(n_cg = `Not treated`, n_tg = Treated)

# contrasts <- contrasts |>
#   mutate(std_err = sqrt((std_cg^2)/n_cg + (std_tg^2)/n_tg))

contrasts <- contrasts |>
  mutate(across(starts_with("std"), ~ .x**2)) |>
  mutate(across(starts_with("std_"),
                ~ .x / get(sub("std_", "n_", cur_column())),
                .names = "ratio_{.col}")) |>
  rowwise() |>
  mutate(std_err = sum(c_across(starts_with("ratio")))) |>
  select(Genus, .data[[config$grouping_variable2]], treatment, difference, starts_with("std"))

to_save[["relabund_for_lm"]] = mO
to_save[["genus_stats"]] = genus_stats
to_save[["diffs"]] = contrasts

if (config$grouping_variable2 != "") { 
  D <- mO %>%
    group_by(Genus,.data[[config$grouping_variable2]]) %>%
    do(tidy(lm(counts ~ .data[[config$grouping_variable1]], data=.))) %>%
    filter(term != "(Intercept)")
} else {
  D <- mO %>%
    group_by(Genus) %>%
    do(tidy(lm(counts ~ .data[[config$grouping_variable1]] + .data[[config$grouping_variable1]], data=.))) %>%
    filter(term != "(Intercept)")
}

D$term <- gsub("treatment","",D$term)
D$term <- gsub("\\.data.*]]","",D$term)
# D$term <- gsub("timepoint","",D$term)

dtbl <- DT::datatable(D, options = list(pageLength=100)) %>% 
  formatStyle('p.value', backgroundColor = styleInterval(config$sig_threshold, c('yellow', 'white')))

## saving HTML file with DT::datatable() output
fname = paste("significant_otus_DT_datatable_", config$suffix, ".html", sep="")
fname = file.path(outdir, "tables", fname)
DT::saveWidget(dtbl, fname)


fname = paste("significant_otus_treatment_within_timepoint_", config$suffix, ".csv", sep="")
fname = file.path(outdir, "tables", fname)
filter(D, p.value <= config$sig_threshold) %>% fwrite(file = fname, sep = ",", col.names = TRUE)

temp <- ungroup(D) |> filter(p.value <= config$sig_threshold)

tmp <- mO |> filter(Genus %in% temp$Genus) |>
  group_by(Genus, .data[[config$grouping_variable1]]) |>
  summarise(avg = mean(counts)) |>
  spread(key = .data[[config$grouping_variable1]], value = avg)

fname = paste("avg_counts_", config$suffix, ".csv", sep="")  
fname = file.path(outdir, "tables", fname)
fwrite(x = tmp, file = fname)

if (config$grouping_variable2 != "") { 
  tmp <- mO |> filter(Genus %in% temp$Genus) |>
    group_by(Genus, .data[[config$grouping_variable2]], .data[[config$grouping_variable1]]) |>
    summarise(avg = mean(counts)) |>
    spread(key = .data[[config$grouping_variable1]], value = avg)
}

if (config$grouping_variable2 != "") {
  
  left_col <- c("Genus", config$grouping_variable2, "treatment")
  right_col <- c("Genus", config$grouping_variable2, "term")
  genus_stats <- genus_stats |>
    inner_join(temp, by = setNames(right_col, left_col))
} else {
  genus_stats <- genus_stats |>
    inner_join(temp, by = c("Genus" = "Genus", "treatment"="term"))
}

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

#################################
### MODEL with time and treatment
#################################
if (config$grouping_variable2 != "") { 
  
  D <- mO %>%
    group_by(Genus) |>
    do(tidy(lm(counts ~ .data[[config$grouping_variable2]] + .data[[config$grouping_variable1]], data=.))) %>%
    filter(term != "(Intercept)")
  
  D$term <- gsub("treatment","",D$term)
  D$term <- gsub("timepoint","",D$term)
  D$term <- gsub("\\.data.*]]","",D$term)
  # D$term <- gsub("timepoint","",D$term)
  
  dtbl <- DT::datatable(D, options = list(pageLength=100)) %>% 
    formatStyle('p.value', backgroundColor = styleInterval(config$sig_threshold, c('yellow', 'white')))
  
  ## saving HTML file with DT::datatable() output
  fname = paste("significant_otus_DT_datatable_time+treat_", config$suffix, ".html", sep="")
  fname = file.path(outdir, "tables", fname)
  DT::saveWidget(dtbl, fname)
  
  
  fname = paste("significant_otus_treatment_across_time_", config$suffix, ".csv", sep="")
  fname = file.path(outdir, "tables", fname)
  filter(D, p.value <= config$sig_threshold) %>% fwrite(file = fname, sep = ",", col.names = TRUE)

}

########################################
## MANUALLY SET THE X VARIABLE IN AOV()
########################################
if (config$grouping_variable2 != "") { 
  
  contrasts <- mO |>
    nest(data = -c(Genus, .data[[config$grouping_variable2]])) |>
    mutate(
      fit = map(data, ~ aov(counts ~ treatment, data = .x)),
      hsd = map(fit, TukeyHSD),
      tidied = map(hsd, tidy)
    ) |>
    unnest(tidied) |>
    select(-c(data,fit,hsd))
} else {
  
  contrasts <- mO |>
    # nest(data = -c(timepoint,metric)) |>
    nest(data = -c(Genus)) |>
    mutate(
      fit = map(data, ~ aov(counts ~ treatment, data = .x)),
      hsd = map(fit, TukeyHSD),
      tidied = map(hsd, tidy)
    ) |>
    unnest(tidied) |>
    select(-c(data,fit,hsd))
}
# vec <- contrasts$Genus %in% temp$Genus
# tmp <- contrasts[vec,]

to_save[["contrasts"]] = contrasts

fname = paste("taxa_contrasts_all_", config$suffix, ".csv", sep="")
fname = file.path(outdir, "tables", fname)
fwrite(contrasts, file = fname)

if (config$grouping_variable2 != "") { 
  
  tmp <- contrasts |>
    filter(!is.na(adj.p.value), adj.p.value < config$sig_threshold) |>
    select(Genus, timepoint, contrast, adj.p.value) |>
    mutate(pvalue = ifelse(adj.p.value <= config$sig_threshold, "<= sig_thr", "> sig_thr")) |>
    select(-adj.p.value)
} else {
  
  tmp <- contrasts |>
    filter(!is.na(adj.p.value), adj.p.value < config$sig_threshold) |>
    select(Genus, contrast, adj.p.value) |>
    mutate(pvalue = ifelse(adj.p.value <= config$sig_threshold, "<= sig_thr", "> sig_thr")) |>
    select(-adj.p.value)
}

colors = c("red","white")

gg <- ggplot(tmp, aes(x = contrast, y = Genus)) + geom_tile(aes(fill=pvalue), color = "white")
if (config$grouping_variable2 != "") gg <- gg + facet_wrap(~timepoint)
gg <- gg + scale_fill_manual(values=colors) + 
  theme(axis.text.x = element_text(angle=90, size = 6),
        axis.text.y = element_text(size = 5),
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

