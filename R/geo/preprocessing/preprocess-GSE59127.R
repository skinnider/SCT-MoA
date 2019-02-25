# Dataset: GSE59127
# Cells: single cells from embryonic day 12.5 kidney
# Values: RPKM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim("data/geo/raw/GSE59127_E12.5_expressed_genes_RPKM.txt.gz")

# fix column names and remove first row
colnames(dat) = dat[1, ]
dat = dat[-1, ]

# remove gene names from the data frame 
genes = dat[, 2]
dat = dat[, 3:88]

# make numeric
dat = as.data.frame(sapply(dat, as.numeric))

# map Entrez to ensembl
dat = map_genes(dat, genes, from = "ENTREZID", to = "ENSEMBL", 
                db = org.Mm.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSE59127.txt")
