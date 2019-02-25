# Dataset: GSE76005
# Cells: embryonic mouse neural progenitor cells
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.csv("data/geo/raw/GSE76005_npc_counts.csv.gz")

# remove gene names from the data frame 
genes = dat[, 1]
dat = dat[, -1]

# map symbol to ensembl
dat = map_genes(dat, genes, from = "SYMBOL", to = "ENSEMBL", db = org.Mm.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSE76005.txt")
