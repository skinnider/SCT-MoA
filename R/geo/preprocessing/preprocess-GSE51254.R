# Dataset: GSE51254
# Cells: single HCT116 cells
# Values: FPKM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
files = list.files("data/geo/raw/GSE51254", full.names = T, pattern = "*.gz")
# analyze only C1 microfluidics cells
files = files[grepl("_C", files)]
accessions = gsub("_.*$", "", basename(files))

# read data
dats = lapply(files, read.delim, header = F)
dats = lapply(dats, function(x) x[, c(1, 10)])
dat = do.call(data.frame, dats)
genes = dat[, 1]

# remove unneeded columns
dat = dat[, grepl("V10", colnames(dat))]
colnames(dat) = accessions

# map symbol to ensembl
dat = map_genes(dat, genes, from = "SYMBOL", to = "ENSEMBL", db = org.Hs.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSE51254.txt")
