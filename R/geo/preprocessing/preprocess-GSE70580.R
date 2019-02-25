# Dataset: GSE70580
# Cells: tonsil innate lymphoid cells
# Values: counts 
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")
library(dplyr)

# read files from ILCs only
files = list.files("data/geo/raw/GSE70580", full.names = T, pattern = "*.gz")
files = files[grepl("_ILC", files)]
accessions = gsub("_.*$", "", basename(files))

# read data
dats = lapply(files, read.delim)
dats = lapply(dats, dplyr::select, X.Gene.symbol, Reads)
dat = do.call(data.frame, dats)

# remove unneeded columns
genes = dat[, 1]
dat = dat[, grepl("Reads", colnames(dat))]

# map to ensembl
dat = map_genes(dat, genes, from = "SYMBOL", to = "ENSEMBL", db = org.Hs.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes
rownames(expr) = accessions

# write
write_and_gzip(expr, "data/geo/processed/GSE70580.txt")
