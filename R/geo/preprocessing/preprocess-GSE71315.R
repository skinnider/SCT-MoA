# Dataset: GSE71315
# Cells: single cells from the developing human neocortex
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim("data/geo/raw/GSE71315_scell_ncounts.genes.thresh.txt.gz")

# analyze single cells only 
genes = dat[, 1]
dat = dat[, grepl("_S", colnames(dat)) | startsWith(colnames(dat), "S")]

# convert to matrix and transpose
expr = t(dat)
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSE71315.txt")
