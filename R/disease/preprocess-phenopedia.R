# Preprocess the Phenopedia database. 
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")
library(tidyverse)
library(org.Hs.eg.db)
library(MeSH.PCR.db)

# load in data
file = "data/disease/phenopedia/Disease-GeneID.txt.gz"
cols = count.fields(file, sep = "\t")
max = max(cols, na.rm = T)
dat = read.delim(file, skip = 3, fill = T, col.names = 1:max)

# remove duplicates, transpose and re-shape
dat = dat[match(unique(dat$X1), dat$X1), ]
dat = as.data.frame(t(dat))
colnames(dat) = dat[1, ]
dat = dat[-1, ]

# clean up to geneIDs, remove blanks
dat = dat %>%
  gather(disease, ENTREZID) %>%
  mutate(ENTREZID = gsub(".*\\(", "", ENTREZID)) %>%
  mutate(ENTREZID = gsub("\\)", "", ENTREZID)) %>%
  filter(ENTREZID != "")

# map to ensembl
map = AnnotationDbi::select(org.Hs.eg.db, keys = unique(dat$ENTREZID),
                     keytypes = "ENTREZID", columns = c("ENTREZID", "ENSEMBL"))
dat = dat %>%
  left_join(map)

# ortholog map to mouse
orthologs = read.delim("data/orthologs/human-mouse-orthologs.txt") %>%
  ## one-to-one
  set_colnames(c("human", "mouse")) %>%
  filter(human != "" & mouse != "") %>%
  group_by(human) %>%
  filter(n_distinct(mouse) == 1) %>%
  ungroup() %>%
  group_by(mouse) %>%
  filter(n_distinct(human) == 1) %>%
  ungroup() %>%
  dplyr::rename(ENSEMBL = human)
dat = dat %>%
  left_join(orthologs, by = 'ENSEMBL') %>%
  dplyr::select(-ENTREZID, -ENSEMBL) %>%
  dplyr::rename(gene = mouse) %>%
  drop_na(gene)

# pull out terms of interest
## cerebrovascular disorders: C0007820
## mental disorders: D001523 (and children thereof)
children = get_children("D001523")
patt = paste0("C0007820|", paste0(children, collapse = "|"))
dat0 = filter(dat, grepl(patt, disease))
write_and_gzip(dat0, "data/disease/phenopedia/phenopedia.txt")
