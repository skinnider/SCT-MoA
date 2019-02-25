# Analyze shared-nearest-neighbor clustering of cell lines with each measure of 
# association.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(scran)
library(igraph)

# create results container
results = data.frame()

# set up parameter values
k_vals = c(2, 5, 10, 20, 50)

# define SNN+louvain clustering function
snn_louvain = function(coexpr, k) {
  # build SNN graph
  ## adapted from scran::buildSNNGraph
  knn = BiocNeighbors::findKNN(-coexpr, k = k, get.distance = F)
  g = .Call(scran:::cxx_build_snn_rank, knn$index)
  edges = g[[1]]
  weights = g[[2]]
  g = make_graph(edges, directed = F)
  E(g)$weight = weights
  g = simplify(g, edge.attr.comb = "first")
  
  # cluster with louvain
  cl = igraph::cluster_louvain(g)
  
  # pull out cluster membership
  clusters = membership(cl)
  return(clusters)
}

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
  
  # analyze a range of k
  for (k in k_vals) {
    # get cluster membership
    clusters = snn_louvain(-coexpr, k)
    # also get ground truth
    types = gsub("_.*$|\\-.*$", "", colnames(coexpr))
    
    # calculate ARI
    ARI = mclust::adjustedRandIndex(clusters, types)
    # calculate NMI
    NMI = ClusterR::external_validation(as.integer(factor(types)), 
                                        as.integer(clusters), method = 'nmi')
    
    # also test within batches
    batches = c("H1", "GM12878")
    patt = paste0(batches, collapse = "|")
    subset = coexpr[grepl(patt, rownames(coexpr)), 
                    grepl(patt, rownames(coexpr))]
    clusters = snn_louvain(-coexpr, k)
    types = gsub("_.*$|\\-.*$", "", colnames(coexpr))
    ARI_batch = mclust::adjustedRandIndex(clusters, types)
    NMI_batch = ClusterR::external_validation(
      as.integer(factor(types)), as.integer(clusters), method = 'nmi')
    
    # add to results
    results %<>% rbind(list(coefficient = coef, k = k, 
                            ARI = ARI, ARI_batch = ARI_batch,
                            NMI = NMI, NMI_batch = NMI_batch))
  }
}

# write
write.csv(results, "data/clustering/sNN-louvain.csv", row.names = F)

# print results
results %>% group_by(coefficient) %>% 
  summarise(ARI = median(ARI)) %>% arrange(-ARI)
results %>% group_by(coefficient) %>% 
  summarise(NMI = median(NMI)) %>% arrange(-NMI)

# plot
results$coefficient %<>% clean_metric()
labels = levels(with(results, reorder(coefficient, -ARI, median)))
p1 = ggplot(results, aes(x = reorder(coefficient, -ARI, median), y = ARI, 
                         fill = coefficient)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(labels = c(
    expression(rho[p]), labels[2], expression(phi[s]), labels[-c(1:3)])) +
  scale_y_continuous("ARI", limits = c(0, 1), 
                     expand = c(0, 0)) +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) +
  clean_theme + 
  theme(axis.title.x = element_blank(), 
        axis.line.y = element_line(color = 'grey50'),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank())
# p1

# plot NMI
labels = levels(with(results, reorder(coefficient, -NMI, median)))
p2 = ggplot(results, aes(x = reorder(coefficient, -NMI, median), 
                         y = NMI, fill = coefficient)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(labels = c(
    labels[1], expression(rho[p]), expression(phi[s]), labels[-c(1:3)])) +
  scale_y_continuous("NMI", limits = c(0, 1), 
                     expand = c(0, 0)) +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = F) +
  clean_theme + 
  theme(axis.title.x = element_blank(), 
        axis.line.y = element_line(color = 'grey50'),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank())
# p2

# save
s10 = plot_grid(p1, p2, labels = letters, label_size = 10)
ggsave("fig/Figure_S10.pdf", s10, width = 17, height = 8, units = 'cm', 
       useDingbats = F)
