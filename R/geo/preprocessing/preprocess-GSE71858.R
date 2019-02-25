# Dataset: GSE71858
# Cells: radial glia
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.csv("data/geo/raw/GSE71858_ES.gene.counts.csv.gz")

# remove gene names from the data frame 
genes = dat[, 1]
dat = dat[, -1]

# map symbol to ensembl
dat = map_genes(dat, genes, from = "SYMBOL", to = "ENSEMBL", db = org.Hs.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSE71858.txt")
