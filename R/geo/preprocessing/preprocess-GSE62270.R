# Dataset: GSE62270
# Cells: single cells isolated from mouse intestinal organoids
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read whole organoids only (no controls)
files = list.files("data/geo/raw/GSE62270", pattern = "*.gz", full.names = T)
files = files[grepl("Whole_Organoid", files)]

# read data
dats = lapply(files, read.delim)
dat = do.call(data.frame, dats)

# map Entrez to ensembl
genes = gsub("__.*$", "", dat[, 1])
dat = dat[, !grepl("GENEID", colnames(dat))]
dat = map_genes(dat, genes, from = "SYMBOL", to = "ENSEMBL", db = org.Mm.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSE62270.txt")
