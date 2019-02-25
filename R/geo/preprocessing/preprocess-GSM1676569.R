# Dataset: GSM1676569
# Cells: single KMB7 cells
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim("data/geo/raw/GSM1676569_KBM755scmRNA.gene.count.csv.gz")

# remove cell names from the data frame 
genes = gsub("__.*$", "", dat[, 1])
dat = dat[ -1]

# map symbol to ensembl
dat = map_genes(dat, genes, from = "SYMBOL", to = "ENSEMBL", db = org.Hs.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSM1676569.txt")
