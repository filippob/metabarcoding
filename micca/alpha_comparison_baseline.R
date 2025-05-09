
library("broom")
library("ggrepel")
library("ggridges")
library("tidyverse")
library("tidymodels")
library("data.table")
library("gghighlight")

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
    prjfolder = "Documents/cremonesi/tamponi_vaginali",
    analysis_folder = "Analysis",
    conf_file = "Config/mapping_file.csv",
    suffix = "vaginal_swabs",
    nfactors = 1, ## n. of design variables (e.g. treatment and timpoint --> nfactors = 2)
    min_tot_n = 15,
    min_sample = 3,
    project = "", ##! use only for subsetting
    treatment_column = "treatment",
    sample_column = "id",
    subject_column = "CAMPIONI",
    grouping_variable2 = "timepoint",
    grouping_variable1 = "treatment",
    base_treatment = 4, ## reference level within timepoint (e.g. control)
    base_timepoint = 1, ## reference level within treatment (e.g. T0)
    force_overwrite = FALSE
  ))
}


HOME <- Sys.getenv("HOME")
repo = file.path(HOME, config$repo)
prjfolder = file.path(HOME, config$prjfolder)
outdir = file.path(prjfolder,config$analysis_folder)

fname = paste("alpha_comparison_baseline.config_",config$suffix,".RData", sep="")
fname = file.path(outdir, fname)
save(config, file = fname)

## treatment levels as in the metadata file
grouping_variable1 = config$grouping_variable1
grouping_variable2 = config$grouping_variable2

## read metadata
metadata = fread(file.path(prjfolder, config$conf_file))
if (config$treatment_column != "treatment" & sum(grepl(pattern = "treatment", names(metadata))) >= 1) {
  
  metadata <- select(metadata, -treatment)
}
metadata <- metadata |> rename('sample-id' = !!config$sample_column, subject = !!config$subject_column, 
                               treatment = !!config$treatment_column)

if (config$project != "") metadata <- filter(metadata, project == !!config$project)

print(head(metadata))

metadata <- mutate(metadata, `sample-id` = as.character(`sample-id`))
metadata <- metadata |> filter(treatment != "", !is.na(treatment))
if("timepoint" %in% names(metadata)) {
  metadata = metadata |> select(c(`sample-id`, subject, treatment, timepoint)) |> mutate(treatment = as.factor(treatment))
} else metadata = metadata |> select(c(`sample-id`, subject, treatment)) |> mutate(treatment = as.factor(treatment))

## read alpha div data
alpha = fread(file.path(prjfolder, config$analysis_folder, "alpha.csv"))
alpha$`sample-id` = gsub('sample.','',alpha$`sample-id`)
alpha$`sample-id` = gsub('\\.','-',alpha$`sample-id`)

alpha = alpha %>% inner_join(metadata, by = c("sample-id" = "sample-id"))

## remove s.e. columns
alpha <- select(alpha, -c(se.chao1,se.ACE))

## reshaping data
if (grouping_variable2 != "") {
  malpha = gather(alpha, key = "metric", value = "value", -c("sample-id", "subject", all_of(grouping_variable1),all_of(grouping_variable2)))
} else malpha = gather(alpha, key = "metric", value = "value", -c("sample-id", "subject", all_of(grouping_variable1)))

# malpha$treatment = factor(malpha$treatment, levels = c("PC","EU","non-EU+","non-EU-"))

alpha_chg <- malpha %>%
  arrange(metric, subject, timepoint) %>%
  group_by(subject) %>%
  mutate(CHG = value - value[1L]) %>%
  ungroup()

ntimepoints = length(unique(alpha_chg$timepoint))
## 
if (ntimepoints <= 2) {
  
  alpha_chg <- alpha_chg |> filter(timepoint == "T2") |>
    select(-timepoint)
}

## alpha div boxplots
p <- ggplot(alpha_chg, aes(x = .data[[grouping_variable1]], y = CHG))
p <- p + geom_boxplot(aes(color = .data[[grouping_variable1]]), alpha=0.7, width = 0.39, size = 1)
p <- p + geom_jitter(aes(color = .data[[grouping_variable1]]), alpha = 0.8, width = 0.1)
p <- p + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "green","darkred"))
if (grouping_variable2 != "") {
  p <- p + facet_wrap(~metric, scales = "free")
} else p <- p + facet_wrap(~metric, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 0))
p

fname = paste("alpha_boxplot_bsl_chg", config$suffix, ".png", sep="")
fname = file.path(outdir, "figures", fname)
ggsave(filename = fname, plot = p, device = "png", width = 7.5, height = 8, dpi = 300)

## density plots
if (ntimepoints <= 2) {
  
  q <- ggplot(alpha_chg, aes(y=.data[[grouping_variable1]], x=CHG, fill=.data[[grouping_variable1]])) +
    geom_density_ridges(scale=0.9, alpha = 0.5) +
    theme(legend.position="none") + facet_wrap(~metric, scales = "free")  
} else {

  q <- ggplot(alpha_chg, aes(y=.data[[grouping_variable2]], x=CHG, fill=.data[[grouping_variable1]])) +
    geom_density_ridges(scale=0.9, alpha = 0.5) +
    theme(legend.position="none") + facet_grid(.data[[grouping_variable2]]~metric, scales = "free")
  # q
} 

print(q)

fname = paste("alpha_density_bsl-chg_", config$suffix, ".png", sep="")
fname = file.path(outdir, "figures", fname)
ggsave(filename = fname, plot = q, device = "png", width = 5, height = 7)

### Linear model
base_treatment = levels(as.factor(alpha$treatment))[config$base_treatment]
base_timepoint = levels(as.factor(alpha$timepoint))[config$base_timepoint]

treats = unique(alpha_chg$treatment)
treats = c(base_treatment, as.character(treats[treats != base_treatment]))
alpha_chg$treatment = factor(alpha_chg$treatment, levels = treats)

if (config$nfactors == 1) {
  
  D <- alpha_chg %>%
    group_by(metric) %>%
    do(tidy(lm(CHG ~ .data[[grouping_variable1]], data = .))) %>%
    filter(term != "(Intercept)", term != "Residuals")
  
  df <- alpha_chg %>%
    group_by(metric) %>%
    do(tidy(anova(lm(CHG ~ .data[[grouping_variable1]], data = .)))) %>%
    filter(term != "(Intercept)", term != "Residuals")
} else {
  
  D <- alpha_chg %>%
    group_by(metric) %>%
    do(tidy(lm(CHG ~ .data[[grouping_variable2]] + .data[[grouping_variable1]], data = .))) %>%
    filter(term != "(Intercept)", term != "Residuals")
  
  df <- alpha_chg %>%
    group_by(metric) %>%
    do(tidy(anova(lm(CHG ~ .data[[grouping_variable2]] + .data[[grouping_variable1]], data = .)))) %>%
    filter(term != "(Intercept)", term != "Residuals")
}

filter(D, p.value < 0.05) |>
  print()

if (ntimepoints > 2) {
  
  alpha_stats <- alpha_chg %>% 
    group_by(metric, .data[[grouping_variable2]], .data[[grouping_variable1]]) %>%
    summarise(avg = round(mean(value),3), std = round(sd(value),3))
} else {
  
  alpha_stats <- alpha_chg %>% 
    group_by(metric, .data[[grouping_variable1]]) %>%
    summarise(avg = round(mean(value),3), std = round(sd(value),3))
}

alpha_stats[alpha_stats$avg > 10,"avg"] <- round(alpha_stats$avg[alpha_stats$avg > 10],1)
alpha_stats[alpha_stats$std > 1,"std"] <- round(alpha_stats$std[alpha_stats$std > 1],1)

dd <- alpha_stats |>
  unite("avg_", avg:std, sep = "+/-") |>
  spread(key = treatment, value = avg_)

fname = paste("alpha__treatment_stats_bsl-chg", config$suffix, ".csv", sep="")
dir.create(file.path(outdir, "tables"), showWarnings = FALSE)
fname = file.path(outdir, "tables", fname)
fwrite(dd, file = fname)

if (ntimepoints > 2) {
  
  avg_timepoint <- alpha_chg |>
    group_by(metric, timepoint) |>
    summarise(avg = mean(CHG)) |>
    rename(term = timepoint)
  
  temp <- avg_timepoint |>
    rename(baseline = avg) |>
    filter(term == base_timepoint) |>
    select(-term)
  
  diff_timepoint <- avg_timepoint |>
    left_join(temp, by = "metric") |>
    mutate(diff = avg - baseline, pct_change = diff/baseline*100)
}

avg_treatment <- alpha_chg |>
  group_by(metric, treatment) |>
  summarise(avg = mean(CHG)) |>
  rename(term = treatment)

temp <- avg_treatment |>
  rename(baseline = avg) |>
  filter(term == base_treatment) |>
  select(-term)

diff_treatment <- avg_treatment |>
  left_join(temp, by = "metric") |>
  mutate(diff = avg - baseline, pct_change = diff/baseline*100)

if (ntimepoints > 2) {
  
  avg <- avg_timepoint |> bind_rows(avg_treatment)
  diff <- diff_timepoint |> bind_rows(diff_treatment) |>
    select(-c(avg,baseline))
} else {
  
  avg <- avg_treatment
  diff <- diff_treatment
}

D$term  <- gsub("\\.data.*\\]","",D$term)

D <- D |>
  left_join(avg, by = c("metric" = "metric", "term" = "term"))

D <- D |>
  left_join(select(diff, -avg), by = c("metric" = "metric", "term" = "term"))


fname = paste("alpha_significance_treatment_bsl-chg", config$suffix, ".csv", sep="")
dir.create(file.path(outdir, "tables"), showWarnings = FALSE)
fname = file.path(outdir, "tables", fname)
fwrite(D, file = fname)

levels(D$metric) <- c("chao","ace","fisher","n_otu","shannon","simpson","equit.","simps_e")

D <- D |> mutate(Color = ifelse(p.value < 0.05, "red", "gray"))

df$term = gsub("\\.data\\[*\"","",df$term)
df$term = gsub("\\]","",df$term)
df$term = gsub("\"","",df$term)

w0 <- ggplot(df, aes(x=term, y=p.value, colour = metric))
w0 <- w0 + geom_point(position=position_dodge(width = 0.5)) 
w0 <- w0 + gghighlight(p.value < 0.05, label_key = metric, label_params = list(size = 4, label.size = 0.5, max.overlaps = 2)) 
w0 <- w0 + geom_hline(yintercept=0.05, linetype="dashed", color = "red", alpha = 0.3)
w0 <- w0 + xlab("term") + coord_cartesian(ylim = c(-0.01,0.1))
w0 <- w0 + theme(axis.text.x = element_text(angle=0, size = 9),
               strip.text.x = element_text(size = 12),
               axis.title = element_text(size = 11))
w0

# D$term <- gsub("Treated","TG",D$term)

w <- ggplot(D, aes(x=term, y=p.value))
w <- w + geom_jitter(aes(colour=metric), size = 2, width = 0.2)
w <- w + gghighlight(p.value < 0.05, label_key = metric, label_params = list(size = 4, label.size = 0.5, max.overlaps = 10))
w <- w + geom_hline(yintercept=0.05, linetype="dashed", color = "red", alpha = 0.3)
w <- w + xlab("term") + coord_cartesian(ylim = c(-0.01,1))
w <- w + theme(axis.text.x = element_text(angle=0, size = 9),
               strip.text.x = element_text(size = 12),
               axis.title = element_text(size = 11))
w


w1 <- ggplot(D, aes(x = term, y = p.value, label = metric))
w1 <- w1 + geom_jitter(aes(group=metric, colour=Color), size = 2, width = 0.2)
w1 <- w1 + geom_hline(yintercept=0.05, linetype="dashed", color = "red", alpha = 0.3)
w1 <- w1 + scale_color_identity() + theme_bw()
w1 <- w1 + xlab("term") + coord_cartesian(ylim = c(-0.01,1))
w1 <- w1 + theme(axis.text.x = element_text(angle=0, size = 9),
                 strip.text.x = element_text(size = 12),
                 axis.title = element_text(size = 11))
w1 <- w1 + geom_label_repel(data = subset(D, p.value < 0.05),
                            size = 3,
                            max.overlaps = Inf,
                            box.padding   = 0.5,
                            point.padding = 0.5,
                            force = 1,
                            segment.size  = 0.2,
                            color = "red")
w1


dir.create(file.path(outdir, "figures"), showWarnings = FALSE)

fname = paste("alpha_significance_bsl-chg_", config$suffix, ".png", sep="")
fname = file.path(outdir, "figures", fname)
ggsave(filename = fname, plot = w1, device = "png", width = 9, height = 7)

########################################
## MANUALLY SET THE X VARIABLE IN AOV()
########################################

if("timepoint" %in% names(alpha_chg)) {
  
  contrasts <- alpha_chg |>
    # nest(data = -c(timepoint,metric)) |>
    nest(data = -c(metric, timepoint)) |>
    mutate(
      fit = map(data, ~ aov(CHG ~ treatment, data = .x)),
      hsd = map(fit, TukeyHSD),
      tidied = map(hsd, tidy)
    ) |>
    unnest(tidied) |>
    select(-c(data,fit,hsd))
} else {
  
  contrasts <- alpha_chg |>
    # nest(data = -c(timepoint,metric)) |>
    nest(data = -c(metric)) |>
    mutate(
      fit = map(data, ~ aov(CHG ~ treatment, data = .x)),
      hsd = map(fit, TukeyHSD),
      tidied = map(hsd, tidy)
    ) |>
    unnest(tidied) |>
    select(-c(data,fit,hsd))
}

fname = paste("alpha_contrasts_bsl-chg_", config$suffix, ".csv", sep="")
fname = file.path(outdir, "tables", fname)
fwrite(contrasts, file = fname)

colors = c("red","white")

if("timepoint" %in% names(alpha_chg)) {
  
  tmp <- contrasts |>
    filter(!is.na(adj.p.value)) |>
    select(metric, timepoint, contrast, adj.p.value) |>
    mutate(pvalue = ifelse(adj.p.value <= 0.05, "<= 0.05", "> 0.05")) |>
    select(-adj.p.value)
  
  gg <- ggplot(tmp, aes(x = contrast, y = metric)) + geom_tile(aes(fill=pvalue), color = "white") + 
    facet_wrap(~timepoint) + 
    scale_fill_manual(values=colors) + 
    theme(axis.text.x = element_text(angle=90, size = 6),
          legend.key = element_rect(color="black"),
          panel.background = element_rect(fill = 'white', color = 'white'),
          panel.grid.major = element_line(color = 'white'),
          panel.grid.minor = element_line(color = 'white'))
} else {
  
  tmp <- contrasts |>
    filter(!is.na(adj.p.value)) |>
    select(metric, contrast, adj.p.value) |>
    mutate(pvalue = ifelse(adj.p.value <= 0.05, "<= 0.05", "> 0.05")) |>
    select(-adj.p.value)
  
  gg <- ggplot(tmp, aes(x = contrast, y = metric)) + geom_tile(aes(fill=pvalue), color = "white") + 
    scale_fill_manual(values=colors) + 
    theme(axis.text.x = element_text(angle=90, size = 6),
          legend.key = element_rect(color="black"),
          panel.background = element_rect(fill = 'white', color = 'white'),
          panel.grid.major = element_line(color = 'white'),
          panel.grid.minor = element_line(color = 'white'))
}


fname = paste("alpha_contrasts_significant_", config$suffix, ".png", sep="")
fname = file.path(outdir, "figures", fname)
ggsave(filename = fname, plot = gg, device = "png")

to_save <- list("alpha_sig_plot_bsl-chg"=w1, "contrasts_bsl-chg_plot"=gg, "alpha_stats_bsl-chg"=dd, "alpha_sign_stats_bsl-chg" = D, "contrasts_tbl_bsl-chg"=contrasts)

fname = paste("alpha_results_bsl-chg_", config$suffix, ".RData", sep="")
fname = file.path(outdir, fname)

save(to_save, file = fname)


