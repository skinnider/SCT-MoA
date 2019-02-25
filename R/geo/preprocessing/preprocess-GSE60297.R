# Dataset: GSE60297
# Cells: mouse medullary thymic epithelial cells
# Values: FPKM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim("data/geo/raw/GSE60297_copy_number.txt.gz")

# remove gene names from the data frame 
genes = dat[, 1]
dat = dat[, -1]

# convert to matrix and transpose
expr = t(dat)
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSE60297.txt")
