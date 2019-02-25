# Plot and statistically test differences in median AUC between measures of 
# association as a function of gene expression filtering.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
source("R/theme.R")
source("R/functions.R")

# read data
dat = read.delim("data/function/auroc_filtered.txt.gz")
# get median for each dataset
med = dat %>%
  group_by(dataset, coefficient, threshold) %>%
  dplyr::summarise(median_auc = median(auroc)) %>%
  ungroup() %>%
  mutate(dataset = gsub("-expr.*$", "", dataset))
# clean up coefficient
med$coefficient %<>% clean_metric() %>% as.character()

# plot, faceted by threshold
p1 = ggplot(med, aes(x = reorder_within(
  coefficient, -median_auc, threshold, median), y = median_auc, 
  fill = coefficient)) + 
  facet_wrap(~ paste0(threshold, "%"), ncol = 3, scales = 'free_x') + 
  geom_hline(aes(yintercept = 0.5), linetype = 'dotted') +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) +
  scale_x_reordered() + 
  scale_y_continuous("AUC", limits = c(0.45, 0.65), expand = c(0, 0)) + 
  theme_sc +
  theme(axis.title.x = element_blank(),
        axis.line.y = element_line(color = 'grey50'),
        panel.grid.major.x = element_line(color = 'grey80', 
                                          linetype = 'dotted'),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = 'grey90', color = 'white'))
p1
ggsave("fig/Figure_S12.pdf", p1, width = 18, height = 16, units = 'cm', 
       useDingbats = F)

# plot, faceted by coefficient
p2 = ggplot(med, aes(x = factor(threshold), y = median_auc, 
                     fill = coefficient)) + 
  facet_wrap(~ coefficient, ncol = 6, scales = 'free_x',
             labeller = label_wrap_gen()) + 
  geom_hline(aes(yintercept = 0.5), linetype = 'dotted') +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) +
  scale_x_discrete("maximum % of dropouts") +
  scale_y_continuous("AUC", limits = c(0.45, 0.65), expand = c(0, 0)) + 
  theme_sc +
  theme(axis.line.y = element_line(color = 'grey50'),
        panel.grid.major.x = element_line(color = 'grey80', 
                                          linetype = 'dotted'),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        strip.background = element_rect(fill = 'grey90', color = 'white'),
        strip.text = element_text(size = 6))
p2
ggsave("fig/Figure_S13.pdf", p2, width = 17, height = 13, units = 'cm', 
       useDingbats = F)
