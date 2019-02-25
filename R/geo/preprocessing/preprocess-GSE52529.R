# Dataset: GSE52529
# Cells: Single differentiating primary human myoblasts
# Values: FPKM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim("data/geo/raw/GSE52529_fpkm_matrix.txt.gz")
# remove version numbers
rownames(dat) = gsub("\\..*$", "", rownames(dat))
# filter out genes with zeroes across the board
dat = dat[rowSums(dat) > 0, ]
# transpose
expr = t(dat)
# write
write_and_gzip(expr, "data/geo/processed/GSE52529.txt")
