# Dataset: GSE81076
# Cells: single cells from human pancreatic islets
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim("data/geo/raw/GSE81076_D2_3_7_10_17.txt.gz")

# map symbols to Ensembl
symbols = gsub("__.*$", "", dat[, 1])
map = AnnotationDbi::select(org.Hs.eg.db, keys = symbols, keytype = 'SYMBOL',
                            columns = 'ENSEMBL')
genes = map$ENSEMBL[match(symbols, map$SYMBOL)]

# convert to matrix and transpose
expr = t(dat[, -1])
colnames(expr) = genes
expr = expr[, !is.na(colnames(expr))]

# write
write_and_gzip(expr, "data/geo/processed/GSE81076.txt")
