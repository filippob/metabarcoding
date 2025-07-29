
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

print("PHYLUM COMPOSITION OF THE MICROBIOME")

## PARAMETERS
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
    prjfolder = "Documents/Suikerbiet/its_2025/its_kyo",
    analysis_folder = "Analysis",
    otu_norm_file = "otu_norm_CSS.csv",
    conf_file = "Config/mapping_file.csv",
    suffix = "its-kyo",
    nfactors = 1, ## n. of design variables (e.g. treatment and timpoint --> nfactors = 2)
    min_tot_n = 10,
    min_sample = 2,
    project = "",
    sample_column = "sample-id",
    sample_prefix = "",
    treatment_column = "Type",
    grouping_variable2 = "timepoint",
    grouping_variable1 = "treatment",
    exp_levels = paste(c("Susceptible", "Resistant"), collapse = ","), ## !! THE FIRST LEVEL IS THE BENCHMARK !! Not treated,Treated
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
if("timepoint" %in% names(metadata)) {
  metadata = metadata |> select(c(`sample-id`, treatment, timepoint)) |> mutate(treatment = as.factor(treatment))
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

if("timepoint" %in% names(metadata)) {
  temp = select(metadata, c(sample, treatment, timepoint))
} else temp = select(metadata, c(sample, treatment, !!config$grouping_variable2))

mO <- mO %>% inner_join(temp, by = "sample")

if("timepoint" %in% names(metadata)) {
  
  D <- mO %>% 
    group_by(treatment,timepoint,Phylum) %>% 
    dplyr::summarise(avg = sum(counts))
  
  D <- mO |> group_by(treatment, timepoint) |>
    mutate(tot = sum(counts)) |>
    group_by(treatment,timepoint, Phylum) |>
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

p <- ggplot(D, aes(x=factor(1), y=avg, fill=Phylum)) + geom_bar(width=1,stat="identity", alpha = 0.8)
p <- p + coord_polar(theta='y')
if (config$grouping_variable2 != "") {
  p <- p + facet_grid(.data[[config$grouping_variable2]]~.data[[config$grouping_variable1]])
} else p <- p + facet_wrap(~treatment)
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
fwrite(x = arrange(D, treatment, desc(avg)), file = fname, sep = ",")

to_save <- list("phylum_relabund"=D)
