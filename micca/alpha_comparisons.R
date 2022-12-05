library("broom")
library("tidyverse")
library("data.table")
library("gghighlight")

alpha = fread("Analysis/microbiome_chlorine/micca/results/alpha.csv")
metadata = fread("Config/mapping_file_chlorine.csv")
alpha$`sample-id` = gsub('sample.','',alpha$`sample-id`)
alpha$`sample-id` = gsub('\\.','-',alpha$`sample-id`)

metadata = select(metadata, c(`sample-id`, treatment, type))
alpha = alpha %>% inner_join(metadata, by = "sample-id")

## remove s.e. columns
alpha <- select(alpha, -c(se.chao1,se.ACE))

malpha = gather(alpha, key = "metric", value = "value", -c("sample-id","treatment","type"))

p <- ggplot(malpha, aes(x = treatment, y = value))
p <- p + geom_boxplot(aes(color = treatment), alpha=1, width = 0.33, size = 1)
p <- p + geom_jitter(aes(color = treatment), alpha = 0.4)
p <- p + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
p <- p + facet_grid(metric~type, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90))
p

ggsave(filename = "Analysis/microbiome_chlorine/micca/results/figures/alpha_boxplot.png", plot = p, device = "png", width = 7, height = 10)

D <- malpha %>%
  group_by(type,metric) %>%
  do(tidy(lm(value ~ treatment, data = .))) %>%
  filter(term != "(Intercept)")

alpha_stats <- malpha %>% 
  group_by(type,metric) %>%
  summarise(avg = mean(value), std = sd(value))

D <- alpha_stats %>% inner_join(D, by = c("type", "metric"))
fwrite(D, file = "Analysis/microbiome_chlorine/micca/results/alpha_significance.csv")

levels(D$metric) <- c("chao","ace","fisher","n_otu","shannon","simpson","equit.","simps_e")

w <- ggplot(D, aes(x=type, y=p.value))
w <- w + geom_jitter(aes(group=metric, colour=metric), width = 0.2)
w <- w + gghighlight(p.value < 0.05, label_params = list(size = 3, label.size = 0.10, max.overlaps = 20))
w <- w + geom_hline(yintercept=0.05, linetype="dashed", color = "red", alpha = 0.2)
w <- w + xlab("type")
w <- w + theme(axis.text.x = element_text(angle=90, size = 7),
               strip.text.x = element_text(size = 12),
               axis.title = element_text(size = 11))
w

ggsave(filename = "Analysis/microbiome_chlorine/micca/results/figures/alpha_significance.png", plot = w, device = "png", width = 9, height = 7)
