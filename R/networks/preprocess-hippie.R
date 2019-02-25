# Preprocess PPI networks from HIPPIE.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/functions.R")

# read original network
ppi = read.delim("data/networks/HIPPIE/raw.txt.gz") %>%
  dplyr::select(protein.A, protein.B)

# map to Ensembl (human)
map = read_tsv("data/identifiers/HUMAN_9606_idmapping.dat.gz",
               col_names = c("uniprot", "db", "id"))
ens = dplyr::filter(map, db == 'Ensembl') %>%
  dplyr::select(-db)
mapped = ppi %>% 
  dplyr::rename(uniprot = protein.A) %>% 
  left_join(ens, by = 'uniprot') %>% 
  dplyr::select(-uniprot) %>%
  dplyr::rename(protein.A = id, uniprot = protein.B) %>%
  left_join(ens, by = 'uniprot') %>% 
  dplyr::select(-uniprot) %>%
  dplyr::rename(protein.B = id) %>%
  drop_na(protein.A, protein.B)
## remove duplicates
mapped %<>% group_by(protein.A, protein.B) %>%
  dplyr::slice(1) %>%
  ungroup()
# write
write_and_gzip(mapped, "data/networks/HIPPIE/human.txt")

# map one-to-one mouse orthologs 
orthologs = read.delim("data/orthologs/human-mouse-orthologs.txt", 
                       col.names = c("human", "mouse")) %>%
  dplyr::filter(!is.na(human) & !is.na(mouse) & human != "" & mouse != "")
## subset to one-to-one orthologs
single_human = names(which(table(orthologs$human) == 1))
single_mouse = names(which(table(orthologs$mouse) == 1))
orthologs %<>% dplyr::filter(human %in% single_human & 
                               mouse %in% single_mouse)
# map PPI network to mouse
mouse = mapped %>% 
  dplyr::rename(human = protein.A) %>% 
  left_join(orthologs, by = 'human') %>% 
  dplyr::select(-human) %>%
  dplyr::rename(protein.A = mouse, human = protein.B) %>%
  left_join(orthologs, by = 'human') %>% 
  dplyr::select(-human) %>%
  dplyr::rename(protein.B = mouse) %>%
  drop_na(protein.A, protein.B)
## remove duplicates (check)
mouse %<>% group_by(protein.A, protein.B) %>%
  dplyr::slice(1) %>%
  ungroup()
# write
write_and_gzip(mouse, "data/networks/HIPPIE/mouse.txt")
