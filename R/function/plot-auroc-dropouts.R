# Plot the relationship between median AUC and % zeroes across measures of 
# association.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")
source("R/functions.R")

# read data
dat = read.delim("data/function/auroc.txt.gz")
# clean coefficients
dat$coefficient %<>% clean_metric()
# fix loom dataset names
dat %<>% mutate(dataset = gsub("-expr.*$", "", dataset))

# get median AUC for each dataset/coefficient
med = dat %>%
  group_by(dataset, coefficient) %>%
  summarise(median_AUC = median(auroc)) %>%
  ungroup()

# add in information about dropouts
prop = read.csv("data/geo/datasets.csv")
med %<>% dplyr::rename(Dataset = dataset) %>%
  left_join(prop, by = 'Dataset')

# filter to one per publication
opp = read.delim("data/geo/one-dataset-per-publication.txt")
med0 = med %>%
  filter(Dataset %in% opp$Dataset)

# plot median AUC vs. % zeroes
pal = colorRampPalette(rev(brewer.pal(n = 12, name = "Paired")))(17)
labels = levels(factor(med$coefficient))
p1 = ggplot(med0, aes(x = Dropouts, y = median_AUC, color = coefficient)) + 
  geom_point(size = 0.8) + 
  geom_smooth(method = 'lm', alpha = 0, size = 0.6) + 
  scale_y_continuous("AUC") + 
  scale_color_manual(values = pal, name = '', labels = c(
    labels[-c(16:17)], expression(rho[p]), expression(phi[s]))) +
  clean_theme + 
  guides(color = guide_legend(ncol = 4, byrow = F)) +
  theme(legend.position = 'bottom', 
        legend.text.align = 0,
        axis.line = element_line(color = 'grey50'))
p1

# test 
rhos = data.frame()
for (coef in unique(med0$coefficient)) {
  subset = dplyr::filter(med0, coefficient == coef)
  cor = cor.test(subset$Dropouts, subset$median_AUC, method = 'spearman')
  rhos %<>% rbind(list(coef = coef, rho = cor$statistic, p = cor$p.value))
}
rhos %>% arrange(p)

# plot median AUC vs. # cells
p2 = ggplot(med0, aes(x = Cells, y = median_AUC, color = coefficient)) + 
  geom_point(size = 0.8) + 
  geom_smooth(method = 'lm', alpha = 0, size = 0.6) + 
  scale_y_continuous("AUC") + 
  scale_color_manual(values = pal, name = '', labels = c(
    labels[-c(16:17)], expression(rho[p]), expression(phi[s]))) +
  scale_x_log10(breaks = c(1, 10, 100, 1e3, 1e4),
                labels = trans_format('log10', math_format(10^.x))) + 
  annotation_logticks(sides = 'b') +
  guides(color = guide_legend(ncol = 4, byrow = F)) +
  clean_theme + 
  theme(legend.position = 'bottom', 
        legend.text.align = 0,
        axis.line = element_line(color = 'grey50'))
p2

# test 
rhos = data.frame()
for (coef in unique(med0$coefficient)) {
  subset = dplyr::filter(med0, coefficient == coef)
  cor = cor.test(subset$Cells, subset$median_AUC, method = 'spearman')
  rhos %<>% rbind(list(coef = coef, rho = cor$statistic, p = cor$p.value))
}
rhos %>% arrange(p) %>% mutate(padj = p.adjust(p, 'BH'), sig = padj < 0.05)

# plot together 
row1 = plot_grid(p1 + theme(legend.position = 'none'),
                 p2 + theme(legend.position = 'none'), 
                 labels = letters, label_size = 10)
row2 = get_legend(p2)
p = plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 0.3))
p
ggsave("fig/Figure_S2.pdf", width = 17.4, height = 11, units = 'cm', 
       useDingbats = F)

# plot by protocol
protocols = c("CEL-Seq2", "Chromium", "Drop-seq", "GemCode", "inDrop",
              "Smart-seq2", "SMARTer", "SMARTer C1", "STRT-seq")
p3 = ggplot(med %>% filter(Protocol %in% protocols),
            aes(x = reorder_within(coefficient, -median_AUC, Protocol, median), 
                y = median_AUC, fill = coefficient)) + 
  scale_x_reordered() + 
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Protocol, scales = 'free_x') + 
  geom_hline(aes(yintercept = 0.5), linetype = 'dotted') +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) +
  scale_y_continuous("AUC", limits = c(0.45, 0.65), expand = c(0, 0)) + 
  theme_sc + 
  theme(axis.title.x = element_blank(),
        axis.line.y = element_line(color = 'grey50'),
        panel.grid.major.x = element_line(color = 'grey80', 
                                          linetype = 'dotted'),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = 'grey90', color = 'white'),
        strip.text = element_text(size = 6))
p3
ggsave("fig/Figure_S3.pdf", p3, width = 17.4, height = 19, units = 'cm',
       useDingbats = F)
## fix labels manually

# plot by coverage
## from: https://teichlab.github.io/scg_lib_structs/
med0 %<>% mutate(coverage = ifelse(
  Protocol %in% c("Smart-seq2", "SMARTer C1", "SMARTer"), 'full-length', ifelse(
    Protocol == 'STRT-seq', '5\'', '3\'')))
p4 = ggplot(filter(med0, Protocol != ""),
            aes(x = reorder_within(coefficient, -median_AUC, coverage, median), 
                y = median_AUC, fill = coefficient)) + 
  scale_x_reordered() + 
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ coverage, scales = 'free_x') + 
  geom_hline(aes(yintercept = 0.5), linetype = 'dotted') +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) +
  scale_y_continuous("AUC", limits = c(0.45, 0.65), expand = c(0, 0)) + 
  theme_sc + 
  theme(axis.title.x = element_blank(),
        axis.line.y = element_line(color = 'grey50'),
        panel.grid.major.x = element_line(color = 'grey80', 
                                          linetype = 'dotted'),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = 'grey90', color = 'white'))
p4
ggsave("fig/Figure_S4.pdf", p4, width = 17.4, height = 8, units = 'cm',
       useDingbats = F)
## fix labels manually
