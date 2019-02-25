# Dataset: GSE70844
# Cells: single cortical neurons sequenced after patch-clamp recording
# Values: FPKM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.csv("data/geo/raw/GSE70844_Fuzik_et_al_molcounts.csv.gz")

# remove gene names from the data frame 
genes = dat[, 1]
dat = dat[, -1]

# map symbol to ensembl
dat = map_genes(dat, genes, from = "ALIAS", to = "ENSEMBL", db = org.Mm.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSE70844.txt")
