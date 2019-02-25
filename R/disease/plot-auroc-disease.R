# Plot and statistically test differences in AUC between measures of
# association derived from Phenopedia disease gene EGAD analysis.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/functions.R")
source("R/theme.R")

# read data
files = list.files("data/disease/auroc", pattern = ".gz", full.names = T)
dat = files %>% 
  map(read.delim) %>% 
  bind_rows() %>%
  mutate(term = gsub(" \\(Psychology)", "-Psychology", term)) %>%
  separate(term, c("name", "CUI"), sep = "\\(") %>%
  mutate(CUI = gsub("\\)", "", CUI))

# filter to CNS disorders in mouse brain atlas
children = get_children("D001523")
psych = filter(dat, CUI %in% children & grepl("-expr", dataset))

# get median for each cell type
med = psych %>%
  group_by(dataset, coefficient) %>%
  dplyr::summarise(median_auc = median(auroc, rm.na = T)) %>%
  ungroup()
# clean up coefficient
med$coefficient %<>% clean_metric() %>% as.character()
# plot
labels1 = levels(with(med, reorder(coefficient, -median_auc, median)))
p1 = ggplot(med, aes(x = stats::reorder(coefficient, -median_auc, median), 
                     y = median_auc, fill = coefficient)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(aes(yintercept = 0.5), linetype = 'dotted') +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) +
  scale_x_discrete(labels = c(expression(rho[p]), expression(phi[s]),
                              labels1[-c(1:2)])) +
  scale_y_continuous("AUC", limits = c(0.4, 0.75), expand = c(0, 0)) +
  theme_sc +
  theme(axis.title.x = element_blank(),
        axis.line.y = element_line(color = 'grey50'),
        panel.grid.major.x = element_line(color = 'grey80',
                                          linetype = 'dotted'),
        axis.ticks.x = element_blank())
p1

# filter to cerebrovascular disorders in brain vasculature atlas
vasc = filter(dat, grepl("GSE98816", dataset) & CUI == "C0007820")

# plot
cvd = vasc %>%
  mutate(dataset = gsub("^.*_", "", dataset),
         coefficient = as.character(clean_metric(coefficient)))
labels2 = setdiff(unique(cvd$coefficient), c("ϕs", "ρp"))
opt = 'lapaz'
p2 = ggplot(cvd, aes(x = reorder(dataset, auroc, median), 
                      y = coefficient, fill = auroc)) + 
  geom_tile() + 
  scale_fill_scico('AUC', palette = opt, limits = c(0.39, 0.6501),
                   breaks = seq(0.45, 0.65, 0.1)) +
  scale_x_discrete("\n", expand = c(0, 0)) +
  scale_y_discrete('', expand = c(0, 0),
                   labels = c(labels2, expression(rho[p]), expression(phi[s]))) +
  coord_flip() +
  theme_sc +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'top',
        legend.direction = 'horizontal',
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.spacing = margin(rep(0, 4)),
        legend.box.margin = margin(rep(0, 4)),
        legend.box.spacing = unit(0, 'cm'))
p3 = ggplot(cvd, aes(x = reorder(dataset, auroc, median), 
                      y = auroc, fill = '1')) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = scico(1, palette = opt, begin = 0.75), guide = F) +
  scale_y_continuous('AUC', expand = c(0, 0)) +
  # coord_cartesian(ylim = c(0.4, 0.65)) +
  coord_flip(ylim = c(0.45, 0.6001)) +
  theme_sc +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(color = 'grey50'),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = 'grey80',
                                          linetype = 'dotted'))
# plot together 
p4 = p2 + p3 + plot_layout(ncol = 2, widths = c(4, 1.25))

# plot separately and combine manually 
ggsave("fig/Figure_5a.pdf", p1, width = 8.9, height = 8, units = 'cm')
ggsave("fig/Figure_5b.pdf", p4, width = 8.9, height = 7.5, units = 'cm')
