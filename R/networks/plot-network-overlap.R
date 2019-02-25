# Plot overlap between different networks and single-cell coexpression networks. 
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(patchwork)
source("R/theme.R")
source("R/functions.R")

# read data
dat = read.delim("data/networks/overlap.txt.gz") %>%
  mutate(dataset = gsub("-expr.*$", "", dataset))

# randomly keep a single network from each publication
keep = read.delim("data/geo/one-dataset-per-publication.txt")[[1]]
dat %<>% filter(dataset %in% keep)

# plot default cutoff (50k edges)
dat1 = filter(dat, cutoff == 5e4) %>%
  mutate(coefficient = clean_metric(coefficient))
labeler = as_labeller(function(x) forcats::fct_recode(
  x, Signalling = "OmniPath", "Metabolic pathways" = "Reactome", 
  "Text mining" = "STRING", "PPIs" = "HIPPIE"))
labels = with(dat1, reorder(coefficient, z_score, median)) %>% levels()
p1 = ggplot(dat1, aes(x = reorder(coefficient, z_score, median), 
                      y = z_score, fill = coefficient)) + 
  facet_wrap(~ network, scales = 'free_x', labeller = labeler, ncol = 4) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous("Z score") +
  scale_x_discrete(labels = c(
    labels[-c(16:17)], expression(rho[p]), expression(phi[s]))) +
  scale_fill_manual(name = '', values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) + 
  coord_flip() + 
  clean_theme + 
  theme(legend.position = 'right',
        axis.text.y = element_text(angle = 0, hjust = 1),
        panel.grid.major.y = element_line(color = 'grey80', 
                                          linetype = 'dotted'),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color = 'grey50'),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = 'grey90', color = 'white'))
p1
ggsave("fig/Figure_2.pdf", p1, width = 18.3, height = 9, units = 'cm')

# plot 20k and 100k edges
dat2 = filter(dat, cutoff == 2e4) %>%
  mutate(coefficient = clean_metric(coefficient))
labels2 = with(dat2, reorder(coefficient, z_score, median)) %>% levels()
p2 = ggplot(dat2, aes(x = reorder(coefficient, z_score, median), 
                      y = z_score, fill = coefficient)) + 
  facet_wrap(~ network, scales = 'free_x', labeller = labeler, ncol = 4) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous("Z score") +
  scale_x_discrete(labels = c(
    labels2[-c(15:17)], expression(rho[p]), expression(phi[s]), labels2[17])) +
  scale_fill_manual(name = '', values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) + 
  coord_flip() + 
  clean_theme + 
  theme(legend.position = 'right',
        axis.text.y = element_text(angle = 0, hjust = 1),
        panel.grid.major.y = element_line(color = 'grey80', 
                                          linetype = 'dotted'),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color = 'grey50'),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = 'grey90', color = 'white'))
p2

dat3 = filter(dat, cutoff == 1e5) %>%
  mutate(coefficient = clean_metric(coefficient))
labels3 = with(dat3, reorder(coefficient, z_score, median)) %>% levels()
p3 = ggplot(dat3, aes(x = reorder(coefficient, z_score, median), 
                      y = z_score, fill = coefficient)) + 
  facet_wrap(~ network, scales = 'free_x', labeller = labeler, ncol = 4) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous("Z score") +
  scale_x_discrete(labels = c(
    labels3[-c(16:17)], expression(rho[p]), expression(phi[s]))) +
  scale_fill_manual(name = '', values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) + 
  coord_flip() + 
  clean_theme + 
  theme(legend.position = 'right',
        axis.text.y = element_text(angle = 0, hjust = 1),
        panel.grid.major.y = element_line(color = 'grey80', 
                                          linetype = 'dotted'),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color = 'grey50'),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = 'grey90', color = 'white'))
p3

# plot together
p4 = plot_grid(p2, p3, labels = letters, ncol = 1, label_size = 10)
ggsave("fig/Figure_S6.pdf", p4, width = 18.3, height = 18, units = 'cm')
