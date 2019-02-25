# Dataset: GSE52583
# Cells: Single lung alveolar epithelial cells
# Values: FPKM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
files = list.files("data/geo/raw/GSE52583", pattern = "*.gz", full.names = T)
# ignore three bulk/no-cell controls
files = files[!grepl("GSM1271882|GSM1271883|GSM1271944", files)]
accessions = gsub("_.*$", "", basename(files))

# read data
dats = lapply(files, read.delim, header = F)
dats = lapply(dats, function(x) x[, c(1, 10)])
dat = do.call(data.frame, dats)
genes = dat[, 1]

# remove unneeded columns
dat = dat[, grepl("V10", colnames(dat))]
colnames(dat) = accessions

# map symbol to ensembl
dat = map_genes(dat, genes, from = "SYMBOL", to = "ENSEMBL", db = org.Mm.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSE52583.txt")
