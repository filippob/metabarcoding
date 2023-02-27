
library("broom")
library("tidyverse")
library("tidymodels")
library("data.table")
library("gghighlight")

## PARAMETERS
HOME <- Sys.getenv("HOME")
repo = file.path(HOME, "Documents/cremonesi/metabarcoding")
prj_folder = file.path(HOME, "Documents/cremonesi/suini_bontempo/pig_feces")
analysis_folder = "Analysis/results"
fname = "filtered_otu/otu_table_filtered.biom"
conf_file = "Config/zinco_poroso_rectum.csv"
outdir = file.path(analysis_folder)
suffix = "zn_porous_feces"

alpha = fread(file.path(prj_folder, analysis_folder, "alpha.csv"))
metadata = fread(file.path(prj_folder,conf_file))
metadata <- metadata |> select(`progressivo MISEQ`, `GRUPPO sperimentale`, TIMEPOINT, Progetto, tipologia)
metadata <- rename(metadata, `sample-id` = `progressivo MISEQ`, treatment = `GRUPPO sperimentale`, timepoint = TIMEPOINT) |> 
  mutate(`sample-id` = as.character(`sample-id`))
alpha$`sample-id` = gsub('sample.','',alpha$`sample-id`)
alpha$`sample-id` = gsub('\\.','-',alpha$`sample-id`)

metadata = select(metadata, c(`sample-id`, timepoint, treatment))
alpha = alpha %>% inner_join(metadata, by = "sample-id")

## remove s.e. columns
alpha <- select(alpha, -c(se.chao1,se.ACE))

malpha = gather(alpha, key = "metric", value = "value", -c("sample-id","timepoint","treatment"))

p <- ggplot(malpha, aes(x = treatment, y = value))
p <- p + geom_boxplot(aes(color = treatment), alpha=0.9, width = 0.33, size = 1)
p <- p + geom_jitter(aes(color = treatment), alpha = 0.3)
p <- p + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "green"))
p <- p + facet_grid(metric~timepoint, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90))
p

fname = paste("alpha_boxplot_", suffix, ".png", sep="")
fname = file.path(prj_folder, outdir, "figures", fname)
ggsave(filename = fname, plot = p, device = "png", width = 7, height = 9)

D <- malpha %>%
  group_by(metric, timepoint) %>%
  do(tidy(lm(value ~ treatment, data = .))) %>%
  filter(term != "(Intercept)")

alpha_stats <- malpha %>% 
  group_by(metric, timepoint, treatment) %>%
  summarise(avg = round(mean(value),3), std = round(sd(value),3))

temp <- alpha_stats |> select(-std) |> spread(key = "treatment", value = avg) |>
  mutate(diffT1 = T1-CTR, diffT2 = T2-CTR, diffT3 = T3-CTR) |>
  select(-c(CTR,T1,T2,T3)) |>
  rename(T1 = diffT1, T2 = diffT2, T3 = diffT3) |>
  gather(key = "treatment", value = "difference_vs_ctrl", -c(metric,timepoint))

alpha_stats <- temp |> 
  inner_join(alpha_stats, by = c("metric" = "metric", "timepoint" = "timepoint", "treatment" = "treatment")) |>
  select(-avg)

D$term <- gsub("treatment","",D$term)

D <- alpha_stats %>% inner_join(D, by = c("metric" = "metric", "timepoint" = "timepoint", "treatment" = "term"))
fname = paste("alpha_significance_", suffix, ".csv", sep="")
fname = file.path(prj_folder, outdir, "tables", fname)
fwrite(D, file = fname)

levels(D$metric) <- c("chao","ace","fisher","n_otu","shannon","simpson","equit.","simps_e")

w <- ggplot(ungroup(D), aes(x=timepoint, y=p.value))
w <- w + geom_jitter(aes(group=metric, colour=metric), width = 0.2)
w <- w + gghighlight(p.value < 0.05, label_params = list(size = 2, label.size = 0.10, max.overlaps = 20))
w <- w + geom_hline(yintercept=0.05, linetype="dashed", color = "red", alpha = 0.2)
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

contrasts <- malpha |>
  nest(data = -c(timepoint,metric)) |>
  mutate(
    fit = map(data, ~ aov(value ~ treatment, data = .x)),
    hsd = map(fit, TukeyHSD),
    tidied = map(hsd, tidy)
  ) |>
  unnest(tidied) |>
  select(-c(data,fit,hsd))

fname = paste("alpha_contrasts_", suffix, ".csv", sep="")
fname = file.path(prj_folder, outdir, "tables", fname)
fwrite(D, file = fname)
