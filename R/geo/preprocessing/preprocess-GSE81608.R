# Dataset: GSE81608
# Cells: single cells from human pancreatic islets
# Values: RPKM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")
library(tidyverse)
library(magrittr)

# read data
dat = read.delim("data/geo/raw/GSE81608/GSE81608_human_islets_rpkm.txt.gz")

# read gene and sample annotations
samples = read.delim("data/geo/raw/GSE81608/human_islet_cell_identity.txt")
genes = read.csv("data/geo/raw/GSE81608/human_gene_annotation.csv")

# match genes to symbols
dat %<>% left_join(genes, by = 'gene.id') %>%
  dplyr::select(symbol, everything())

# map symbols to Ensembl
symbols = dat$symbol
map = AnnotationDbi::select(org.Hs.eg.db, keys = symbols, keytype = 'SYMBOL',
                            columns = 'ENSEMBL')
genes = map$ENSEMBL[match(symbols, map$SYMBOL)]

# convert to matrix and transpose
expr = t(dat[, -c(1:2)])
colnames(expr) = genes
expr = expr[, !is.na(colnames(expr)) & !duplicated(colnames(expr))]

# split up by cell types
cell_types = unique(samples$cell.type)
cell_types = cell_types[!grepl("contaminated", cell_types)]
for (type in cell_types) {
  cells = gsub("\\..*$", "", samples$raw.file[samples$cell.type == type])
  subset = expr[rownames(expr) %in% cells,]
  output_file = paste0("data/geo/processed/GSE81608_", type, ".txt")
  write_and_gzip(subset, output_file)
  message("wrote cell type `", type, "` with ", nrow(subset), " cells")
}
