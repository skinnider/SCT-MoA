# Dataset: GSE67980
# Cells: single prostate cancer circulating tumor cells 
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim("data/geo/raw/GSE67980_readCounts.txt.gz")

# subset to single CTCs
genes = as.character(dat[, 2])
dat = dat[, 6:127]

# map to ensembl
dat = map_genes(dat, genes, from = "ENTREZID", to = "ENSEMBL", 
                db = org.Hs.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSE67980.txt")
