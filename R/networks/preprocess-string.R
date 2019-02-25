# Preprocess text-mining gene-gene associations from STRING. 
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(rtracklayer)

# read human
hs = read.table("data/networks/STRING/9606.protein.links.detailed.v10.5.txt.gz", 
                skip = 1)
hs %<>% dplyr::select(V1, V2, V9) %>%
  dplyr::rename(id1 = V1, id2 = V2, textmining = V9) %>%
  dplyr::filter(textmining >= 500) %>%
  dplyr::mutate(id1 = gsub("9606\\.", "", id1),
                id2 = gsub("9606\\.", "", id2))

# map Ensembl protein to genes
gtf = as.data.frame(import("data/ensembl/Homo_sapiens.GRCh38.91.gtf.gz"))
g2p = gtf %>% dplyr::filter(type == 'CDS') %>%
  dplyr::select(gene_id, protein_id) %>%
  group_by(gene_id, protein_id) %>%
  dplyr::slice(1) %>%
  ungroup()
map = with(g2p, setNames(gene_id, protein_id))
hs$id1 = map[hs$id1]
hs$id2 = map[hs$id2]

# remove duplicates
hs %<>% group_by(id1, id2) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(-textmining)

# write 
write_and_gzip(hs, "data/networks/STRING/human.txt")

# read mouse
mm = read.table("data/networks/STRING/10090.protein.links.detailed.v10.5.txt.gz",
                skip = 1)
mm %<>% dplyr::select(V1, V2, V9) %>%
  dplyr::rename(id1 = V1, id2 = V2, textmining = V9) %>%
  dplyr::filter(textmining >= 500) %>%
  dplyr::mutate(id1 = gsub("10090\\.", "", id1),
                id2 = gsub("10090\\.", "", id2))

# map Ensembl protein to genes
gtf = as.data.frame(import("data/ensembl/Mus_musculus.GRCm38.91.gtf.gz"))
g2p = gtf %>% dplyr::filter(type == 'CDS') %>%
  dplyr::select(gene_id, protein_id) %>%
  group_by(gene_id, protein_id) %>%
  dplyr::slice(1) %>%
  ungroup()
map = with(g2p, setNames(gene_id, protein_id))
mm$id1 = map[mm$id1]
mm$id2 = map[mm$id2]

# remove duplicates
mm %<>% group_by(id1, id2) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(-textmining)

# write 
write_and_gzip(mm, "data/networks/STRING/mouse.txt")
