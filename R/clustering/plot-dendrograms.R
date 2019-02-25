# Plot cell line dendrograms from hierarchical clustering with each measure of 
# association.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(ggtree)
library(ape)
source("R/theme.R")

# create list of dendrograms
dendrograms = list()

# analyze each matrix in turn 
files = list.files("data/clustering/matrices", pattern = "*.Rdata", 
                   full.names = T)
for (file in files) {
  coef = gsub("^.*-|\\.Rdata", "", basename(file))
  load(file, verbose = T) ## coexpr
  
  # replace missing values with median (ZI kendall)
  coexpr[is.na(coexpr)] = median(coexpr, na.rm = T)
  # replace infinite values with the minimum (binomial)
  coexpr[is.infinite(coexpr)] = min(coexpr, na.rm = T)
  
  # scale metrics to [-1, 1] range
  if (!dismay::is_bounded(coef)) {
    coexpr = scales::rescale(coexpr, to = c(-1, 1))
  }
  
  # hierarchically cluster
  clust = hclust(as.dist(-coexpr))
  # get cell type labels
  labs = clust$labels[clust$order]
  labels1 = gsub("_.*$|\\-.*$|\\..*$", "", labs)
  # get batch labels
  labels2 = ifelse(grepl("_B", labs), 
                   gsub("^.*_|-.*$|\\..*$", "", labs), NA)

  # plot tree
  p1 = ggtree(as.phylo(clust), branch.length = "none", size = 0.2) +
    scale_x_reverse() + 
    ggtitle(ifelse(coef == 'phi_s', expression(phi[s]), ifelse(
      coef == 'rho_p', expression(rho[p]), clean_metric(coef)))) +
    coord_flip() +
    theme(plot.margin = margin(rep(0, 4)),
          plot.title = element_text(size = 6, hjust = 0.5))
  p1
  
  # plot heatmap
  hm = data.frame(x = clust$labels, idx = seq_along(clust$labels),
                  cell_type = labels1, batch = labels2)
  long = gather(hm, 'y', 'fill', 3:4)
  # first row: cell type  
  pal1 = scico::scico(palette = 'batlow', n = 7)
  p2 = ggplot(hm, aes(x = idx, y = '1', fill = cell_type)) + 
    geom_tile() +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = pal1, name = 'cell type', guide =  F) +
    theme_sc + 
    theme(axis.text = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'right',
          legend.direction = 'horizontal',
          plot.margin = margin(rep(0, 4)))
  # p2
  # second row: batch
  pal2 = c(scico::scico(palette = 'roma', n = 10)[c(3, 8)])
  p3 = ggplot(hm, aes(x = idx, y = '1', fill = batch)) + 
    geom_tile() + 
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = pal2, name = 'batch', guide = F,
                      na.value = 'grey90') +
    theme_sc + 
    theme(axis.text = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'right',
          legend.direction = 'horizontal',
          plot.margin = margin(rep(0, 4)))
  # p3
  
  # plot together
  p = p1 + p2 + p3 + plot_layout(ncol = 1, heights = c(8, 1, 1))
  # p
  # save
  ggsave(paste0("fig/panels/dendrogram-", which(files == file), ".pdf"),
         p, width = 4, height = 3, units = 'cm', useDingbats = F)
  
  # add to list
  dendrograms[[coef]] = p
}

# finally, plot legend
l1 = p2 + 
  scale_fill_manual(values = pal1, name = 'cell type') +
  theme(legend.direction = 'horizontal',
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6))
l2 = p3 + 
  scale_fill_manual(values = pal2, name = 'batch', na.value = 'grey90') +
  theme(legend.direction = 'horizontal',
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6))
legend = l1 + l2 + plot_layout(ncol = 1)
legend
ggsave("fig/panels/dendrogram-legend.pdf", legend, width = 10, height = 10, 
       units = 'cm', useDingbats = F)

## assemble manually - patchwork can't handle
