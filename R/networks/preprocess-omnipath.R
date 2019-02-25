# Preprocess the OmniPath signalling pathway network.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# read data
net = read.delim("data/networks/OmniPath/interactions.txt.gz")

# map IDs to human
map = read_tsv("data/identifiers/HUMAN_9606_idmapping.dat.gz",
               col_names = c("uniprot", "db", "id"))
ens = dplyr::filter(map, db == 'Ensembl') %>%
  dplyr::select(-db)
mapped = net %>% 
  dplyr::rename(uniprot = source) %>% 
  left_join(ens, by = 'uniprot') %>% 
  dplyr::select(-uniprot) %>%
  dplyr::rename(id1 = id, uniprot = target) %>%
  left_join(ens, by = 'uniprot') %>% 
  dplyr::select(-uniprot) %>%
  dplyr::rename(id2 = id) %>%
  drop_na(id1, id2)
## remove duplicates
mapped %<>% group_by(id1, id2) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(id1, id2)
# write
write_and_gzip(mapped, "data/networks/OmniPath/human.txt")

# map one-to-one mouse orthologs 
orthologs = read.delim("data/orthologs/human-mouse-orthologs.txt", 
                       col.names = c("human", "mouse")) %>%
  dplyr::filter(!is.na(human) & !is.na(mouse) & human != "" & mouse != "")
## subset to one-to-one orthologs
single_human = names(which(table(orthologs$human) == 1))
single_mouse = names(which(table(orthologs$mouse) == 1))
orthologs %<>% dplyr::filter(human %in% single_human & 
                               mouse %in% single_mouse)
# map signalling network to mouse
mouse = mapped %>% 
  dplyr::rename(human = id1) %>% 
  left_join(orthologs, by = 'human') %>% 
  dplyr::select(-human) %>%
  dplyr::rename(id1 = mouse, human = id2) %>%
  left_join(orthologs, by = 'human') %>% 
  dplyr::select(-human) %>%
  dplyr::rename(id2 = mouse) %>%
  drop_na(id1, id2)
## remove duplicates (check)
mouse %<>% group_by(id1, id2) %>%
  dplyr::slice(1) %>%
  ungroup()
# write
write_and_gzip(mouse, "data/networks/OmniPath/mouse.txt")
