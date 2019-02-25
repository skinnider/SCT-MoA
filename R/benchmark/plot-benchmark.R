# Plot wall clock time to construct coexpression networks with each measure
# of association.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
source("R/functions.R")
source("R/theme.R")

# read data
dat = list.files("data/benchmark", pattern = "*.txt", full.names = T) %>%
  map(read.delim) %>%
  bind_rows() %>% 
  dplyr::rename(median_time = median) %>%
  mutate(median_time = median_time / 10^9)
# clean up coefficients
dat$coef %<>% clean_metric() %>% as.character()
# get labels 
labels = levels(with(dat, reorder(coef, median_time, median)))
# plot 
breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1e3)
p = ggplot(dat, aes(x = reorder(coef, median_time, median), y = median_time,
                    fill = coef)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) +
  scale_x_discrete(labels = c(
    labels[1:6], expression(rho[p]), expression(phi[s]), 
    labels[9:17])) +
  scale_y_log10("Time (s)", expand = c(0, 0),
                breaks = breaks, labels = as.character(breaks)) +
  annotation_logticks(sides = 'l', short = unit(0.05, 'cm'), 
                      mid = unit(0.1, 'cm'), long = unit(0.15, 'cm')) +
  theme_sc + 
  theme(axis.title.x = element_blank(),
        axis.line.y = element_line(color = 'grey50'),
        panel.grid.major.x = element_line(color = 'grey80', 
                                          linetype = 'dotted'),
        axis.ticks.x = element_blank())
p
ggsave("fig/Figure_S11.pdf", p, width = 8.9, height = 8, units = 'cm')
