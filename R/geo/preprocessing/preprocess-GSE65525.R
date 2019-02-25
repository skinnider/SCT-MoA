# Dataset: GSE65525
# Cells: single K562 cells
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.csv("data/geo/raw/GSE65525/GSM1599500_K562_cells.csv.bz2")

# map symbols to Ensembl genes
genes = dat[, 1]
dat = map_genes(dat[, -1], genes, from = "SYMBOL", to = "ENSEMBL",
                db = org.Hs.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSE65525.txt")
