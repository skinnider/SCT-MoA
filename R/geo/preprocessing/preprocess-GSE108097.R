# Dataset: GSE108097
# Cells: single cells from mouse tissues
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/functions.R")

# read cell type assignments
## downloaded from https://figshare.com/s/865e694ad06d5857db4b
types = read.csv("data/geo/raw/GSE108097/MCA_cellassignments.csv.gz") %>%
  dplyr::filter(grepl('brain', Tissue, ignore.case = T))

# read cells
files = list.files("data/geo/raw/GSE108097", pattern = "*_Brain*", 
                   full.names = T)
dats = map(files, ~ read.table(., check.names = F))
dat = map(dats, ~ rownames_to_column(., 'gene')) %>%
  Reduce(function(x, y) full_join(x, y, by = 'gene'), .)

# map symbols to Ensembl
symbols = dat$gene
map = AnnotationDbi::select(org.Mm.eg.db, keys = symbols, keytype = 'SYMBOL',
                            columns = 'ENSEMBL')
genes = map$ENSEMBL[match(symbols, map$SYMBOL)]

# convert to matrix
expr = t(dat[, -1])
colnames(expr) = genes
expr = expr[, !is.na(colnames(expr)) & !duplicated(colnames(expr))]

# write each cell type separately 
expr = expr[rownames(expr) %in% types$Cell.name, ]
types %<>% dplyr::filter(Cell.name %in% rownames(expr))
minor_types = names(which(table(types$ClusterID) >= 50))
for (type in minor_types) {
  cells = types$Cell.name[types$ClusterID == type]
  subset = expr[rownames(expr) %in% cells,]
  output = paste0("data/geo/processed/GSE108097_", type, ".txt")
  write_and_gzip(subset, output)
  message("wrote cell type ", type, " with ", nrow(subset), " cells")
}
