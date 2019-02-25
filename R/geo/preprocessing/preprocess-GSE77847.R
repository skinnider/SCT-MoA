# Dataset: GSE77847
# Cells: single cells from a mouse model of acute myeloid leukemia
# Values: TPM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim("data/geo/raw/GSE77847_Mouse-scRNASeq-AML_RSEM-TPM.txt.gz")

# remove gene names from the data frame 
genes = dat[, 1]
dat = dat[, -1]

# map symbol to ensembl
dat = map_genes(dat, genes, from = "SYMBOL", to = "ENSEMBL", db = org.Mm.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSE77847.txt")
