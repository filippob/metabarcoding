## libraries
library("dplyr")
library("tools")
library("vegan")
library("plotly")
library("phyloseq")
library("metagenomeSeq")

## PARAMETERS
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 1) {
  
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
    repo = "Documents/cremonesi/metabarcoding",
    prjfolder = "Documents/cremonesi/vitelli_giulia_sala_2025",
    analysis_folder = "Analysis",
    conf_file = "Config/mapping_file.csv",
    suffix = "calves_colostrum",
    nfactors = 2, ## n. of design variables (e.g. treatment and timpoint --> nfactors = 2)
    min_tot_n = 15,
    min_sample = 3,
    project = "",
    treatment_column = "treatment",
    factor_cov = "timepoint",
    sample_id = "", ## from metadata
    force_overwrite = FALSE
  ))
}

## SET UP
HOME <- Sys.getenv("HOME")
repo = file.path(HOME, config$repo)
prjfolder = file.path(HOME, config$prjfolder)
outdir = file.path(prjfolder,config$analysis_folder)

fname = paste("beta-3d.config_",config$suffix,".RData", sep="")
fname = file.path(outdir, fname)
save(config, file = fname)

source(file.path(repo, "support_functions/dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(repo, "support_functions/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/

## READ DATA
## loading data previously imported in phyloseq
fname = file.path(outdir, "phyloseq.RData")
load(fname)

## filter
min_tot_n = config$min_tot_n
min_sample = config$min_sample

# temp = filter_taxa(otu_tax_sample, function(x) sum(x > 1) > (0.20*length(x)), TRUE)
filtered_taxa = filter_taxa(otu_tax_sample, function(x) sum(x > 1) > min_sample, TRUE)
filtered_taxa = filter_taxa(filtered_taxa, function(x) sum(x) > min_tot_n, TRUE)
summary(rowSums(otu_table(filtered_taxa))) |> print()

rm(otu_tax_sample)
gc()

## Preprocessing: e.g. filtering
## making normalization folder
if(!file.exists(file.path(outdir))) dir.create(file.path(outdir), showWarnings = FALSE)

writeLines(" - CSS normalization")
otu_tax_sample_norm = phyloseq_transform_css(filtered_taxa, norm = TRUE, log = FALSE)
sample_data(otu_tax_sample_norm) |> head() |> print()

rm(filtered_taxa)
gc()

## subset for subproject
if(nchar(config$project) > 0) {
  
  otu_norm_subset <- subset_samples(otu_tax_sample_norm, project == config$project)
} else otu_norm_subset <- otu_tax_sample_norm

## changing treatment column
x <- sample_data(otu_tax_sample_norm)
sample_data(otu_norm_subset)$treatment <- dplyr::pull(x, !!config$treatment_column )
sample_data(otu_norm_subset) |> head() |> print()

otus <- as.matrix(otu_table(otu_norm_subset))
otus <- t(otus)
matx <- data.matrix(otus)

### MULTIDIMENSIONAL SCALING
omics.mds = metaMDS(matx, k=3)
omics.scores <- as.data.frame(scores(omics.mds, display = 'site'))

## METADATA AND OMICS PREP
metadata <- sample_data(otu_norm_subset)

if (config$sample_id == "") {
    metadata$id = row.names(metadata)
    metadata <- as_tibble(metadata)
    metadata <- select(metadata, c(id, !!config$treatment_column, !!config$factor_cov))

    omics.scores <- omics.scores |> 
    mutate(id := row.names(omics.scores)) |>
    relocate(id)

    omics.scores <- omics.scores |> inner_join(metadata, by = "id")
    
} else {
    metadata <- as_tibble(metadata)
    metadata <- select(metadata, c(!!config$sample_id, !!config$treatment_column, !!config$factor_cov))

    omics.scores <- omics.scores |> 
    mutate(!!config$sample_i := row.names(omics.scores)) |>
    relocate(!!config$sample_id)

    omics.scores <- omics.scores |> inner_join(metadata, by = !!config$sample_id)
}


## 3D PLOT
temp <- dplyr::select(omics.scores, c(id, NMDS1, NMDS2, NMDS3, !!config$treatment_column))
temp <- temp |> rename(group = !!config$treatment_column)

p <- plot_ly(data = temp, 
        x = ~NMDS1, y = ~NMDS2, z = ~NMDS3,
        type = "scatter3d",
        color = ~ group,
        colors = c('#BF382A', '#0C4B8E', '#B59410')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'NMDS1'),
                      yaxis = list(title = 'NMDS2'),
                      zaxis = list(title = 'NMDS3')),
         annotations = list(
           x = 0.005,
           y = 0.01,
           text = 'Omics clustering',
           xref = 'paper',
           yref = 'paper',
           showarrow = FALSE
         ))


fname = paste("3d-clustr-",config$suffix,"-",config$treatment_column, ".html", sep="")
fname = file.path(outdir, fname)
htmlwidgets::saveWidget(p, fname)

## OPTIONAL second 3D PLOT
if(config$factor_cov != "") {
    
    temp <- dplyr::select(omics.scores, c(id, NMDS1, NMDS2, NMDS3, !!config$factor_cov))
    temp <- temp |> rename(group = !!config$factor_cov)

    p <- plot_ly(data = temp, 
        x = ~NMDS1, y = ~NMDS2, z = ~NMDS3,
        type = "scatter3d",
        color = ~ group,
        colors = c('#BF382A', '#0C4B8E', '#B59410')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'NMDS1'),
                      yaxis = list(title = 'NMDS2'),
                      zaxis = list(title = 'NMDS3')),
         annotations = list(
           x = 0.005,
           y = 0.01,
           text = 'Omics clustering',
           xref = 'paper',
           yref = 'paper',
           showarrow = FALSE
         ))

    fname = paste("3d-clustr-",config$suffix,"-",config$factor_cov, ".html", sep="")
    fname = file.path(outdir, fname)
    htmlwidgets::saveWidget(p, fname)
}

print("DONE!!")