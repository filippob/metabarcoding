
library("broom")
library("tidyverse")
library("data.table")
library("gghighlight")

## PARAMETERS
HOME <- Sys.getenv("HOME")
repo = file.path(HOME, "Documents/cremonesi/metabarcoding")
prj_folder = file.path(HOME, "Documents/cremonesi/suini_bontempo/pig_feces")
analysis_folder = "Analysis/results"
fname = "filtered_otu/otu_table_filtered.biom"
conf_file = "Config/rectum_mapping.csv"
outdir = file.path(analysis_folder)

alpha = fread(file.path(prj_folder, analysis_folder, "alpha.csv"))
metadata = fread(file.path(prj_folder,conf_file))
metadata <- rename(metadata, `sample-id` = id) |> mutate(`sample-id` = as.character(`sample-id`))
alpha$`sample-id` = gsub('sample.','',alpha$`sample-id`)
alpha$`sample-id` = gsub('\\.','-',alpha$`sample-id`)

metadata = select(metadata, c(`sample-id`, timepoint, treatment))
alpha = alpha %>% inner_join(metadata, by = "sample-id")

## remove s.e. columns
alpha <- select(alpha, -c(se.chao1,se.ACE))

malpha = gather(alpha, key = "metric", value = "value", -c("sample-id","timepoint","treatment"))

p <- ggplot(malpha, aes(x = treatment, y = value))
p <- p + geom_boxplot(aes(color = treatment), alpha=1, width = 0.33, size = 1)
p <- p + geom_jitter(aes(color = treatment), alpha = 0.4)
p <- p + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "green"))
p <- p + facet_grid(metric~timepoint, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90))
p

fname = file.path(prj_folder, outdir, "figures/alpha_boxplot.png")
ggsave(filename = fname, plot = p, device = "png", width = 7, height = 9)

D <- malpha %>%
  group_by(metric, timepoint) %>%
  do(tidy(lm(value ~ treatment, data = .))) %>%
  filter(term != "(Intercept)")

alpha_stats <- malpha %>% 
  group_by(metric, timepoint) %>%
  summarise(avg = mean(value), std = sd(value))

D <- alpha_stats %>% inner_join(D, by = c("metric", "timepoint"))
fname = file.path(prj_folder, outdir, "alpha_significance.csv")
fwrite(D, file = fname)

levels(D$metric) <- c("chao","ace","fisher","n_otu","shannon","simpson","equit.","simps_e")

w <- ggplot(D, aes(x=timepoint, y=p.value))
w <- w + geom_jitter(aes(group=metric, colour=metric), width = 0.2)
w <- w + gghighlight(p.value < 0.05, label_params = list(size = 2, label.size = 0.10, max.overlaps = 20))
w <- w + geom_hline(yintercept=0.05, linetype="dashed", color = "red", alpha = 0.2)
w <- w + xlab("term")
w <- w + theme(axis.text.x = element_text(angle=90, size = 7),
               strip.text.x = element_text(size = 12),
               axis.title = element_text(size = 11))
w

fname = file.path(prj_folder, outdir, "figures/alpha_significance.png")
ggsave(filename = fname, plot = w, device = "png", width = 9, height = 7)

