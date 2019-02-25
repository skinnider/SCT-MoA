# Dataset: GSE64960
# Cells: single differentiating cells from the mouse gonad
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.csv("data/geo/raw/GSE64960_single_cell_counts.csv.gz")

# remove bulk sample
dat = dat[, -2]

# remove gene names from the data frame 
genes = dat[, 1]
dat = dat[, -1]

# convert to matrix and transpose
expr = t(dat)
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSE64960.txt")
