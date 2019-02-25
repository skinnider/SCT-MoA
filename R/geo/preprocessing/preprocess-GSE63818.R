# Dataset: GSE63818
# Cells: single primordial germ cells
# Values: FPKM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# analyze PGCs only
files = list.files("data/geo/raw/GSE63818", pattern = "*.gz", full.names = T)
files = files[grepl("PGC", files)] 
## ignore pool/split controls
files = files[!grepl("_ps", files)]

# read data
dats = lapply(files, read.delim)
dat = do.call(data.frame, dats)

# remove unneeded columns
genes = dat[, 1]
genes = Hmisc::capitalize(tolower(genes))
dat = dat[, !grepl("tracking_id", colnames(dat))]

# map to ensembl
dat = map_genes(dat, genes, from = "SYMBOL", to = "ENSEMBL", db = org.Mm.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes

# write
write_and_gzip(expr, "data/geo/processed/GSE63818.txt")
