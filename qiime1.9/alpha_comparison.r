library("broom")
library("tidyverse")
library("data.table")
library("gghighlight")

### PREDIPPING - HERD COMPARISON
alpha = fread("Analysis/microbiome_predipping//qiime1.9/results/alpha.csv")
metadata = fread("Config/mapping_file.csv")
alpha$`sample-id` = gsub('sample.','',alpha$`sample-id`)
alpha$`sample-id` = gsub('\\.','-',alpha$`sample-id`)

metadata = select(metadata, c(`#sampleID`, azienda, treatment, timepoint, herd))
metadata <- rename(metadata, `sample-id` = `#sampleID`) %>% mutate(`sample-id` = as.character(`sample-id`))
alpha$`sample-id` = paste(alpha$`sample-id`,alpha$`sample-id`, sep="")

alpha = alpha %>% inner_join(metadata, by = "sample-id")

## remove s.e. columns
alpha <- select(alpha, -c(se.chao1,se.ACE))

malpha = gather(alpha, key = "metric", value = "value", -c("sample-id","treatment","timepoint","azienda","herd"))

D <- malpha %>% group_by(herd,azienda,treatment,metric) %>%
  summarise(N=n(), avg = round(mean(value),2)) %>%
  spread(key = azienda, value = avg)

fname = "Analysis/microbiome_predipping/qiime1.9/results/alpha_herd.csv"
fwrite(x = D, file = fname, sep = ",")

p <- ggplot(malpha, aes(x = treatment, y = value))
p <- p + geom_boxplot(aes(color = treatment), alpha=1, width = 0.33, size = 1)
p <- p + geom_jitter(aes(color = treatment), alpha = 0.4)
p <- p + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
p <- p + facet_grid(metric~azienda, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, size = 7),
               axis.text.y = element_text(size = 7))
p

fname = "Analysis/microbiome_predipping/qiime1.9/results/alpha_herd.png"
ggsave(filename = fname, plot = p, device = "png", width = 9, height = 8)
