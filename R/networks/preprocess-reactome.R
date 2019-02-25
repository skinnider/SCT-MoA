# Preprocess metabolic pathway co-membership networks from Reactome.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(flavin)

# read Reactome pathways
react = read_tsv("data/networks/reactome/Ensembl2Reactome.txt.gz",
                 col_names = c("ensembl", "pathway", "url", "name", "ec", "sp"))
# filter mouse and human separately
mouse = react %>% dplyr::filter(sp == 'Mus musculus') %>% 
  dplyr::select(ensembl, pathway)
human = react %>% dplyr::filter(sp == 'Homo sapiens') %>% 
  dplyr::select(ensembl, pathway)
# convert to annotation lists
annMm = as_annotation_list(mouse, "ensembl", "pathway") 
annHs = as_annotation_list(human, "ensembl", "pathway")
annMm = annMm[lengths(annMm) >= 2]
annHs = annHs[lengths(annHs) >= 2]
# get all combinations within each list 
pairwise = function(vec) as.data.frame(t(combn(vec, 2)))
mm = map(annMm, pairwise) %>% bind_rows()
hs = map(annHs, pairwise) %>% bind_rows()
colnames(mm) = colnames(hs) = c("gene1", "gene2")
# remove any duplicates
mm %<>% group_by(gene1, gene2) %>% dplyr::slice(1) %>% ungroup()
hs %<>% group_by(gene1, gene2) %>% dplyr::slice(1) %>% ungroup()
# write 
write_and_gzip(mm, "data/networks/reactome/mouse.txt")
write_and_gzip(hs, "data/networks/reactome/human.txt")
