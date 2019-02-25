# Plot the properties of the single-cell coexpression datasets analyzed 
# in this study.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/theme.R")
library(tidyverse)
library(magrittr)

# read manually annotated table (with protocols)
prop = read.csv("data/geo/datasets.csv")

# recategorize very limited protocols
ltd = names(which(table(prop$Protocol) < 5))
prop %<>% mutate(Protocol = ifelse(Protocol %in% ltd, 'other', Protocol))

# recode factor
protocols = setdiff(unique(prop$Protocol), "other")
prop %<>% mutate(Protocol = factor(Protocol, levels = c(protocols, 'other')))

# panel A: number of cells
p1 = ggplot(prop, aes(x = Cells, fill = '1', color = '1')) + 
  geom_histogram(bins = 40, alpha = 0.75) +
  scale_fill_manual(values = c('1' = '#cab2d6'), guide = F) +
  scale_color_manual(values = c('1' = '#cab2d6'), guide = F) +
  scale_x_log10(breaks = c(1, 10, 100, 1e3, 1e4),
                labels = scales::trans_format('log10',math_format(10^.x))) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20), 
                     breaks = seq(0, 20, 5)) +
  labs(x = "Cells", y = "Datasets") + 
  annotation_logticks(sides = 'b', short = unit(0.05, 'cm'),
                      mid = unit(0.1, 'cm'), long = unit(0.15, 'cm'),
                      size = 0.25) +
  theme_sc  + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.line = element_line(color = 'grey50'))
p1

# panel B: number of genes
p2 = ggplot(prop, aes(x = Genes, fill = '1', color = '1')) + 
  geom_histogram(bins = 40, alpha = 0.75) +
  scale_fill_manual(values = c('1' = '#a6cee3'), guide = F) +
  scale_color_manual(values = c('1' = '#a6cee3'), guide = F) +
  scale_x_continuous("Genes (thousands)", labels = function(x) x / 1e3) + 
  scale_y_continuous("Datasets", breaks = seq(0, 10, 2), limits = c(0, 10),
                     expand = c(0, 0)) +
  theme_sc + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.line = element_line(color = 'grey50'))
p2

# panel C: proportion of dropouts
p3 = ggplot(prop, aes(x = Dropouts, fill = '1', color = '1')) + 
  geom_histogram(bins = 40, alpha = 0.75) +
  scale_x_continuous("% zeroes", labels = function(x) x * 100) +
  scale_y_continuous("Datasets", expand = c(0, 0), limits = c(0, 25),
                     breaks = seq(0, 25, 5)) +
  scale_fill_manual(values = c('1' = '#33a02c'), guide = F) +
  scale_color_manual(values = c('1' = '#33a02c'), guide = F) +
  theme_sc + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.line = element_line(color = 'grey50'))
p3

# panel D: dropouts vs. cells 
p4 = ggplot(filter(prop, Protocol != ""), 
            aes(x = Cells, y = Dropouts)) + 
  geom_point(aes(color = Protocol), size = 0.8) + 
  geom_smooth(method = 'lm', alpha = 0.1, linetype = 'dotted', 
              color = 'grey50') +
  scale_color_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(15), name = '') +
  scale_x_log10(breaks = c(1, 10, 100, 1e3, 1e4),
                labels = scales::trans_format('log10',math_format(10^.x))) + 
  
  scale_y_continuous("% zeroes", labels = function(x) 100 * x) +
  annotation_logticks(sides = 'b', short = unit(0.05, 'cm'),
                      mid = unit(0.1, 'cm'), long = unit(0.15, 'cm'),
                      size = 0.25) +
  guides(color = guide_legend(ncol = 1, byrow = F)) +
  theme_sc +
  theme(legend.position = 'right', axis.line = element_line(color = 'grey50'),
        axis.text.x = element_text(angle = 0, hjust = 0.5))
p4

# panel E: sequencing protocols 
p5 = ggplot(prop %>% filter(Protocol != ""),
            aes(x = Protocol, color = '1', fill = '1')) + 
  geom_histogram(stat = 'count', alpha = 0.75) +
  scale_fill_manual(values = c('1' = '#e31a1c'), guide = F) +
  scale_color_manual(values = c('1' = '#e31a1c'), guide = F) +
  scale_y_continuous("Datasets", expand = c(0, 0), limits = c(0, 50)) +
  theme_sc + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(color = 'grey50'),
        axis.title.x = element_blank())
p5

save(p1, p2, p3, p4, p5, file = "fig/Rdata/datasets.Rdata")
