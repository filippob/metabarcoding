

## SET UP
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
    prjfolder = "Documents/Suikerbiet/its_2025",
    analysis_folder = "Analysis",
    otu_norm_file = "otu_norm_CSS.csv",
    conf_file = "Config/mapping_file.csv",
    suffix = "its1f-its2",
    nfactors = 2, ## n. of design variables (e.g. treatment and timpoint --> nfactors = 2)
    min_tot_n = 15,
    min_sample = 3,
    project = "",
    sample_column = "sample-id",
    sample_prefix = "",
    treatment_column = "Type",
    grouping_variable2 = "timepoint",
    grouping_variable1 = "treatment",
    sig_threshold = 0.01,
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

## read metadata
metadata = fread(file.path(prjfolder, config$conf_file))
# if(config$treatment_column != "treatment" & "treatment" %in% names(metadata)) metadata <- rename(metadata, secondary_treatment = treatment)
if(config$treatment_column != "treatment" & "treatment" %in% names(metadata)) metadata$treatment <- NULL
metadata <- metadata |> rename(`sample-id` = !!config$sample_column, treatment = !!config$treatment_column)

if (config$project != "") metadata <- filter(metadata, project == !!config$project)

metadata <- mutate(metadata, `sample-id` = as.character(`sample-id`))
if("timepoint" %in% names(metadata)) {
  metadata = metadata |> select(c(`sample-id`, treatment, timepoint)) |> mutate(treatment = as.factor(treatment))
} else metadata = metadata |> select(c(`sample-id`, treatment, !!config$grouping_variable2))

## convert to factors
metadata <- metadata |> mutate_at(vars(-("sample-id")),as.factor)

## read OTUs (normalised)
otus = fread(otu_norm_fname)
otus <- filter(otus, Phylum != "")


## F:B ratio
names(otus)

## subset otu table (only samples present in the mapping file)
vec <- c(TRUE,vec,rep(TRUE,7))
otus <- as.data.frame(otus)
otus <- otus[,vec]

temp <- otus |>
  select(-c(tax_id,Kingdom, Class, Order, Family, Genus, Species)) |>
  gather(key = "sample", value = "counts", -c("Phylum")) |>
  filter(Phylum != "")

temp <- temp %>% inner_join(select(metadata, -`sample-id`), by = c("sample" = "sample"))

if (config$grouping_variable2 != "") {
  fb <- temp |>
    filter(Phylum %in% c("Firmicutes", "Bacteroidetes")) |>
    group_by(sample, Phylum) |>
    summarise(tot = sum(counts), treatment = unique(treatment), timepoint = unique(timepoint)) |>
    spread(key = Phylum, value = tot) |>
    mutate(FB = Firmicutes/Bacteroidetes, treatment = factor(treatment, levels = exp_levels))
} else {
  
  fb <- temp |>
    filter(Phylum %in% c("Firmicutes", "Bacteroidetes")) |>
    group_by(sample, Phylum) |>
    summarise(tot = sum(counts), treatment = unique(treatment)) |>
    spread(key = Phylum, value = tot) |>
    mutate(FB = Firmicutes/Bacteroidetes, treatment = factor(treatment, levels = exp_levels))
}

to_save[["data_for_fb"]] = temp
to_save[["fb_res"]] = fb

p1 <- ggplot(filter(fb, FB != Inf), aes(x = FB, y = treatment, fill = treatment)) 
p1 <- p1 + geom_density_ridges(aes(point_color = treatment, point_fill = treatment, point_shape = treatment),
                               alpha = .2, point_alpha = 1, jittered_points = TRUE) 
if (config$grouping_variable2 != "") p1 <- p1 + facet_wrap(~timepoint, scales = "free")
p1 <- p1 + scale_point_color_hue(l = 40)
p1 <- p1 + scale_discrete_manual(aesthetics = "point_shape", values = (seq(1,length(exp_levels))+20))
p1

fname = paste("fb_density_", config$suffix, ".png", sep="")
fname = file.path(outdir, fname)
ggsave(filename = fname, plot = p, device = "png")

if (config$grouping_variable2 != "") {
  temp <- filter(fb, FB != Inf) |>
    group_by(timepoint) |>
    do(tidy(lm(FB ~ treatment, data = .))) |>
    filter(term != "(Intercept)") |>
    mutate(term = gsub('treatment','',term)) |>
    rename(treatment = term) |>
    select(timepoint, treatment, p.value)
} else {
  
  temp <- filter(fb, FB != Inf) |>
    do(tidy(lm(FB ~ treatment, data = .))) |>
    filter(term != "(Intercept)") |>
    mutate(term = gsub('treatment','',term)) |>
    rename(treatment = term) |>
    select(treatment, p.value)
}

writeLines(" - FB ratio")
print(temp)

## save results to R object
fname = paste("fb_ratio_results_", config$suffix, ".RData", sep="")
fname = file.path(outdir, fname)
save(to_save, file = fname)

library("DONE!!")
