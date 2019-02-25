# Dataset: GSE69761
# Cells: single cells from the developing mouse lung at E16.5
# Values: FPKM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim("data/geo/raw/GSE69761_36188g-148c-fpkm.txt.gz")

# remove gene names from the data frame 
genes = gsub("\\..*$", "", dat[, 1])
dat = dat[, -1]

# convert to matrix and transpose
expr = t(dat)
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSE69761.txt")
