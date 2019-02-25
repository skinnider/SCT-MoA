# Dataset: GSE74596
# Cells: single mouse natural killer T cells
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
files = list.files("data/geo/raw/GSE74596", full.names = T, pattern = "*.gz")
accessions = gsub("_.*$", "", basename(files))
dats = lapply(files, read.delim, header = F)
dat = do.call(data.frame, dats)

# remove unneeded columns
genes = dat[, 1]
dat = dat[, grepl("V2", colnames(dat))]

# map to ensembl
dat = map_genes(dat, genes, from = "ALIAS", to = "ENSEMBL", db = org.Mm.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes
rownames(expr) = accessions

# write
write_and_gzip(expr, "data/geo/processed/GSE74596.txt")
