# Dataset: GSE71485
# Cells: single hippocampal quiescent neural stem cells 
# Values: TPM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim("data/geo/raw/GSE71485_Single_TPM.txt.gz")

# subset to CFPnuc+/- cells
dat = dat[, startsWith(colnames(dat), "N") | startsWith(colnames(dat), "C")]

# map to Ensembl
dat = map_genes(dat, rownames(dat), from = "ALIAS", to = "ENSEMBL",
                db = org.Mm.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSE71485.txt")
