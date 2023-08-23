
library("broom")
library("tidyverse")
library("tidymodels")
library("data.table")
library("gghighlight")

## PARAMETERS
HOME <- Sys.getenv("HOME")
repo = file.path(HOME, "Documents/cremonesi/metabarcoding")
prj_folder = file.path(HOME, "Documents/cremonesi/suini_bontempo")
analysis_folder = "Analysis/micca/results_zinc_caecum"
# fname = "filtered_otu/otu_table_filtered.biom"
# conf_file = "Config/mapping_file.csv"
conf_file = "Config/zinco_poroso_caecum.csv"
outdir = file.path(prj_folder,analysis_folder)
suffix = "zinc_caecum"

grouping_variable2 = "timepoint"
grouping_variable1 = "treatment"

alpha = fread(file.path(prj_folder, analysis_folder, "alpha.csv"))
metadata = fread(file.path(prj_folder,conf_file))

metadata <- metadata |> select(`progressivo MISEQ`, `GRUPPO sperimentale`, TIMEPOINT, Progetto, tipologia)
# metadata <- metadata |> filter(Progetto == "Zinco Poroso")
metadata <- rename(metadata, `sample-id` = `progressivo MISEQ`, treatment = `GRUPPO sperimentale`, timepoint = TIMEPOINT) |>
  mutate(`sample-id` = as.character(`sample-id`))
metadata = select(metadata, c(`sample-id`, treatment, timepoint)) |> filter(timepoint == "D28")

alpha$`sample-id` = gsub('sample.','',alpha$`sample-id`)
alpha$`sample-id` = gsub('\\.','-',alpha$`sample-id`)

# metadata <- mutate(metadata, id = as.character(id))
# metadata = select(metadata, c(`sample-id`, treatment))
# metadata <- metadata |> select(-c(seq_paola, sampling_date, Tipologia)) |> mutate(id = as.character(id))
# metadata <- metadata |> select(-c(pig)) |> mutate(id = as.character(id))
alpha = alpha %>% inner_join(metadata, by = c("sample-id" = "sample-id"))

## remove s.e. columns
alpha <- select(alpha, -c(se.chao1,se.ACE))

# malpha = gather(alpha, key = "metric", value = "value", -c("sample-id","timepoint","treatment"))
malpha = gather(alpha, key = "metric", value = "value", -c("sample-id",all_of(grouping_variable1),all_of(grouping_variable2)))

p <- ggplot(malpha, aes(x = .data[[grouping_variable1]], y = value))
# p <- p + geom_boxplot(aes(color = substrate), alpha=0.9, width = 0.33, size = 1)
p <- p + geom_jitter(aes(color = .data[[grouping_variable1]]), alpha = 0.8, width = 0.1)
p <- p + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "green","darkred"))
p <- p + facet_grid(metric~.data[[grouping_variable2]], scales = "free")
# p <- p + facet_wrap(~metric, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90))
p

library("ggridges")
q <- ggplot(malpha, aes(y=.data[[grouping_variable1]], x=value, fill=.data[[grouping_variable1]])) +
  geom_density_ridges(scale=0.9, alpha = 0.5) +
  theme(legend.position="none") + facet_wrap(~metric, scales = "free")
q

fname = paste("alpha_plot_", suffix, ".png", sep="")
fname = file.path(prj_folder, outdir, "figures", fname)
ggsave(filename = fname, plot = p, device = "png", width = 5, height = 7)

D <- malpha %>%
  group_by(metric, .data[[grouping_variable2]]) %>%
  # do(tidy(lm(value ~ kit + substrate + treatment, data = .))) %>%
  do(tidy(lm(value ~ .data[[grouping_variable1]], data = .))) %>%
  filter(term != "(Intercept)")

alpha_stats <- malpha %>% 
  group_by(metric, .data[[grouping_variable1]]) %>%
  summarise(avg = round(mean(value),3), std = round(sd(value),3))

# temp <- alpha_stats |> select(-std) |> spread(key = .data[[grouping_variable2]], value = avg) |>
#   mutate(diffA1 = A1-AX, diffA2 = A2-AX) |>
#   select(-c(A1,A2,AX))
# temp <- alpha_stats |> select(-std) |> spread(key = .data[[grouping_variable1]], value = avg) |>
#   mutate(diff = INSETTI-SOIA) |>
#   select(-c(INSETTI,SOIA))
temp <- alpha_stats |> select(-std) |> spread(key = .data[[grouping_variable1]], value = avg) |>
  mutate(diff1 = T1-CTR, diff2 = T2-CTR, diff3 = T3-CTR) |>
  select(-c(CTR,T1,T2,T3))

# 
alpha_stats <- temp |>
  inner_join(alpha_stats, by = c("metric" = "metric"))

# temp <- alpha_stats |> select(-std) |> spread(key = "treatment", value = avg) |>
#   mutate(diffT1 = T1-CTR, diffT2 = T2-CTR, diffT3 = T3-CTR) |>
#   select(-c(CTR,T1,T2,T3)) |>
#   rename(T1 = diffT1, T2 = diffT2, T3 = diffT3) |>
#   gather(key = "treatment", value = "difference_vs_ctrl", -c(metric))

# alpha_stats <- temp |>
#   inner_join(alpha_stats, by = c("metric" = "metric", "kit" = "kit", "substrate" = "substrate")) |>
#   # inner_join(alpha_stats, by = c("metric" = "metric", "treatment" = "treatment")) |>
#   select(-std)

D$term  <- gsub("\\.data.*\\]","",D$term)

# D <- alpha_stats %>% inner_join(D, by = c("metric" = "metric", "timepoint" = "timepoint", "treatment" = "term"))
D <- alpha_stats %>% inner_join(D, by = c("metric" = "metric", setNames("term", grouping_variable1))) |>
  select(-c(avg,std))

fname = paste("alpha_significance_treatment_", suffix, ".csv", sep="")
dir.create(file.path(outdir, "tables"), showWarnings = FALSE)
fname = file.path(outdir, "tables", fname)
fwrite(D, file = fname)

levels(D$metric) <- c("chao","ace","fisher","n_otu","shannon","simpson","equit.","simps_e")

w <- ggplot(ungroup(D), aes(x=.data[[grouping_variable2]], y=p.value))
w <- w + geom_jitter(aes(group=metric, colour=metric), size = 3, width = 0.2)
w <- w + gghighlight(p.value < 0.10, label_params = list(size = 3, label.size = 0.25, max.overlaps = 10))
w <- w + geom_hline(yintercept=0.05, linetype="dashed", color = "red", alpha = 0.3)
w <- w + xlab("term")
w <- w + theme(axis.text.x = element_text(angle=90, size = 7),
               strip.text.x = element_text(size = 12),
               axis.title = element_text(size = 11))
w

fname = paste("alpha_significance_", suffix, ".png", sep="")
fname = file.path(prj_folder, outdir, "figures", fname)
ggsave(filename = fname, plot = w, device = "png", width = 9, height = 7)


## Tukey HSD (contrasts)
writeLines(" - making pairwise contrasts between treatments")
# lmfit = aov(value ~ treatment, data = filter(malpha,metric == "Chao1" & timepoint == "D28"))
# xx = TukeyHSD(lmfit)
# xx$treatment


########################################
## MANUALLY SET THE X VARIABLE IN AOV()
########################################
contrasts <- malpha |>
  # nest(data = -c(timepoint,metric)) |>
  nest(data = -c(metric, timepoint)) |>
  mutate(
    fit = map(data, ~ aov(value ~ treatment, data = .x)),
    hsd = map(fit, TukeyHSD),
    tidied = map(hsd, tidy)
  ) |>
  unnest(tidied) |>
  select(-c(data,fit,hsd))

fname = paste("alpha_contrasts_", suffix, ".csv", sep="")
fname = file.path(outdir, "tables", fname)
fwrite(contrasts, file = fname)

tmp <- contrasts |>
  filter(!is.na(adj.p.value)) |>
  select(metric, timepoint, contrast, adj.p.value) |>
  mutate(pvalue = ifelse(adj.p.value <= 0.05, "<= 0.05", "> 0.05")) |>
  select(-adj.p.value)

colors = c("red","white")

gg <- ggplot(tmp, aes(x = contrast, y = metric)) + geom_tile(aes(fill=pvalue), color = "white") + 
  facet_wrap(~timepoint) + 
  # scale_fill_distiller(palette = "YlGnBu") + 
  scale_fill_manual(values=colors) + 
  theme(axis.text.x = element_text(angle=90, size = 6),
        # axis.text.y = element_text(size = 4),
        # legend.key.size = unit(0.1, 'cm'),
        # strip.text = element_text(size=7),
        # legend.title = element_text(size=6),
        # legend.text = element_text(size=6),
        legend.key = element_rect(color="black"),
        panel.background = element_rect(fill = 'white', color = 'white'),
        panel.grid.major = element_line(color = 'white'),
        panel.grid.minor = element_line(color = 'white'))

fname = paste("alpha_contrasts_significant_", suffix, ".png", sep="")
fname = file.path(outdir, "figures", fname)
ggsave(filename = fname, plot = gg, device = "png")



