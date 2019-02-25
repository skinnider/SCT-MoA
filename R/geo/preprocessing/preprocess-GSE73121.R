# Dataset: GSE73121
# Cells: renal cell carcinoma (primary, metastatic, PDX)
# Values: TPM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")
library(dplyr)

# read files
files = list.files("data/geo/raw/GSE73121", full.names = T, pattern = "*.gz")
# get single-cell files only
files = files[grepl("_SC_", files)]
accessions = gsub("_.*$", "", basename(files))

# read data
dats = lapply(files, read.delim)
dats = lapply(dats, dplyr::select, gene_id, TPM)
dat = do.call(data.frame, dats)

# remove unneeded columns
genes = gsub("\\..*$", "", dat[, 1])
dat = dat[, grepl("TPM", colnames(dat))]

# convert to matrix and transpose
expr = t(dat)
colnames(expr) = genes
rownames(expr) = accessions

# write
write_and_gzip(expr, "data/geo/processed/GSE73121.txt")
