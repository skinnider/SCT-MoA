# Analyze hierarchical clustering of cell lines with each measure of 
# association.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# create results container
results = data.frame()

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
  clusters = cutree(clust, k = 7)
  types = gsub("_.*$|\\-.*$", "", names(clusters))
  
  # calculate ARI
  ARI = mclust::adjustedRandIndex(clusters, types)
  # calculate NMI
  NMI = ClusterR::external_validation(as.integer(factor(types)), 
                                      clusters, method = 'nmi')
  
  # also test within batches
  batches = c("H1", "GM12878")
  patt = paste0(batches, collapse = "|")
  subset = coexpr[grepl(patt, rownames(coexpr)), grepl(patt, rownames(coexpr))]
  clust = hclust(as.dist(-subset))
  clusters = cutree(clust, k = 2)
  types = gsub("_.*$|\\-.*$", "", names(clusters)[grepl(patt, names(clusters))])
  # calculate ARI
  ARI_batch = mclust::adjustedRandIndex(clusters, types)
  # calculate NMI
  NMI_batch = ClusterR::external_validation(as.integer(factor(types)), 
                                            clusters, method = 'nmi')
  
  # add to results
  results %<>% rbind(list(coefficient = coef, ARI = ARI, ARI_batch = ARI_batch,
                          NMI = NMI, NMI_batch = NMI_batch))
}

# write
write.csv(results, "data/clustering/hclust.csv", row.names = F)

# print results
results %>% arrange(-ARI)
results %>% arrange(-ARI_batch)
results %>% arrange(-NMI)
results %>% arrange(-NMI_batch)

# plot main figure 
results$coefficient %<>% clean_metric()
labels = levels(with(results, reorder(coefficient, -ARI)))
p = ggplot(results, aes(x = reorder(coefficient, -ARI), y = ARI, 
                        fill = coefficient)) +
  geom_col() +
  scale_x_discrete(labels = c(expression(rho[p]), expression(phi[s]),
                              labels[-c(1:2)])) +
  scale_y_continuous("adjusted Rand index", limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) +
  clean_theme + 
  theme(axis.title.x = element_blank(), 
        axis.line.y = element_line(color = 'grey50'),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank())
p
ggsave("fig/Figure_4.pdf", p, width = 8.9, height = 8, units = 'cm', 
       useDingbats = F)

# plot within-batch analyses
labels = levels(with(results, reorder(coefficient, -ARI_batch)))
p = ggplot(results, aes(x = reorder(coefficient, -ARI_batch), y = ARI_batch,
                        fill = coefficient)) +
  geom_col() +
  scale_x_discrete(labels = c(labels[1], expression(rho[p]), expression(phi[s]),
                              labels[-c(1:3)])) +
  scale_y_continuous("\nARI, H1/GM12878 cells only", limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) +
  clean_theme + 
  theme(axis.title.x = element_blank(), 
        axis.line.y = element_line(color = 'grey50'),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank())
p
ggsave("fig/Figure_S8.pdf", p, width = 8.9, height = 8, units = 'cm', 
       useDingbats = F)

# plot NMI figure 
labels = levels(with(results, reorder(coefficient, -NMI)))
pA = ggplot(results, aes(x = reorder(coefficient, -NMI), y = NMI, 
                         fill = coefficient)) +
  geom_col() +
  scale_x_discrete(labels = c(expression(rho[p]), expression(phi[s]),
                              labels[-c(1:2)])) +
  scale_y_continuous("NMI", limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) +
  clean_theme + 
  theme(axis.title.x = element_blank(), 
        axis.line.y = element_line(color = 'grey50'),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank())
pA
# within-batch
labels = levels(with(results, reorder(coefficient, -NMI_batch)))
pB = ggplot(results, aes(x = reorder(coefficient, -NMI_batch), y = NMI_batch,
                         fill = coefficient)) +
  geom_col() +
  scale_x_discrete(labels = c(labels[1], expression(rho[p]), expression(phi[s]),
                              labels[-c(1:3)])) +
  scale_y_continuous("NMI", limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) +
  clean_theme + 
  theme(axis.title.x = element_blank(), 
        axis.line.y = element_line(color = 'grey50'),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank())
pB
s9 = plot_grid(pA, pB, labels = letters, label_size = 10)
ggsave("fig/Figure_S9.pdf", s9, width = 17, height = 8, units = 'cm', 
       useDingbats = F)
