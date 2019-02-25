# Plot and statistically test differences in median AUC between measures of 
# association.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
source("R/theme.R")
source("R/functions.R")

# read data
dat = read.delim("data/function/auroc.txt.gz")
# get median for each dataset
med = dat %>%
  group_by(dataset, coefficient) %>%
  dplyr::summarise(median_auc = median(auroc)) %>%
  ungroup() %>%
  mutate(dataset = gsub("-expr.*$", "", dataset))
# clean up coefficeint
med$coefficient %<>% clean_metric() %>% as.character()

# filter to one dataset per publication
opp = read.delim("data/geo/one-dataset-per-publication.txt")
med0 = filter(med, dataset %in% opp$Dataset)

# load other Fig 1 panels
load("fig/Rdata/datasets.Rdata", verbose = T)

# plot
labels = levels(with(med0, reorder(coefficient, -median_auc, median)))
p6 = ggplot(med0, aes(x = reorder(coefficient, -median_auc, median), 
                     y = median_auc, fill = coefficient)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_hline(aes(yintercept = 0.5), linetype = 'dotted') +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) +
  scale_x_discrete(labels = c(expression(rho[p]), expression(phi[s]),
                   labels[-c(1:2)])) +
  scale_y_continuous("AUC", limits = c(0.45, 0.65), expand = c(0, 0)) + 
  theme_sc + 
  theme(axis.title.x = element_blank(),
        axis.line.y = element_line(color = 'grey50'),
        panel.grid.major.x = element_line(color = 'grey80', 
                                          linetype = 'dotted'),
        axis.ticks.x = element_blank(),
        plot.background = element_blank())
p6

# plot together
p4 = p4 +
  guides(color = guide_legend(ncol = 4)) +
  theme(legend.box.margin = margin(rep(0, 4)),
                legend.box.spacing = unit(0, "cm"),
                legend.position = 'bottom')
row1 = plot_grid(p1, p2, p3, p5, ncol = 4, labels = letters,
                 rel_widths = c(1, 1, 1, 1.3), label_size = 10)
row2 = plot_grid(p4, p6, ncol = 2, labels = letters[5:6], 
                 rel_widths = c(1, 1.3), label_size = 10)
fig1 = plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1.9),
                 align = 'center')
ggsave("fig/Figure_1.pdf", fig1, width = 15, height = 11, units = 'cm')

# also plot the entire dataset, without downsampling by publication
labels2 = levels(with(med, reorder(coefficient, -median_auc, median)))
p2 = ggplot(med, aes(x = reorder(coefficient, -median_auc, median),
                     y = median_auc, fill = coefficient)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_hline(aes(yintercept = 0.5), linetype = 'dotted') +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) +
  scale_x_discrete(labels = c(expression(phi[s]), expression(rho[p]),
                              labels2[-c(1:2)])) +
  scale_y_continuous("AUC", limits = c(0.45, 0.65), expand = c(0, 0)) + 
  theme_sc + 
  theme(axis.title.x = element_blank(),
        axis.line.y = element_line(color = 'grey50'),
        panel.grid.major.x = element_line(color = 'grey80', 
                                          linetype = 'dotted'),
        axis.ticks.x = element_blank())
p2
ggsave("fig/Figure_S1.pdf", p2, width = 8.9, height = 8, units = 'cm')

# print median AUCs
med0 %>% 
  group_by(coefficient) %>%
  summarise(median = median(median_auc)) %>%
  arrange(desc(median))

# run statistical tests 
coefs = dismay::metrics()
pvals = matrix(NA, nrow = length(coefs), ncol = length(coefs), 
               dimnames = list(coefs, coefs))
for (i in seq_len(length(coefs))) {
  coef1 = coefs[i]
  message("analyzing coefficient ", coef1, " ...")
  for (j in seq_len(length(coefs))) {
    coef2 = coefs[j]
    if (coef2 == coef1)
      next
    message("  analyzing coefficient ", coef2, " ...")
    # run Brunner--Munzel tests within each dataset
    p = numeric(0)
    for (d in unique(med0$dataset)) {
      x = dat %>% dplyr::filter(dataset == d & coefficient == coef1) %>% 
        pull(auroc)
      y = dat %>% dplyr::filter(dataset == d & coefficient == coef2) %>% 
        pull(auroc)
      test = lawstat::brunner.munzel.test(x, y)
      p[d] = test$p.value
    }
    # Fisher integration 
    fisher = aggregation::fisher(p)
    pvals[coef1, coef2] = fisher
  }
}
write.table(pvals, "data/results/auroc_pvals.txt", quote = F, row.names = T,
            sep = "\t")
