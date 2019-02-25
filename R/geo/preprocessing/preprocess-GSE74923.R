# Dataset: GSE74923
# Cells: activated murine CD8+ T-cell and lymphocytic leukemia cell lines
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim("data/geo/raw/GSE74923_L1210_CD8_processed_data.txt.gz")

# remove gene names
genes = dat[, 1]
dat = dat[, -1]

# map symbol to ensembl
dat = map_genes(dat, genes, from = "SYMBOL", to = "ENSEMBL", db = org.Mm.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes

# write L1210 and CD8 separately
types = c("L1210", "CD8")
for (type in types) {
  subset = expr[grepl(type, rownames(expr)),]
  write_and_gzip(subset, paste0("data/geo/processed/GSE74923_", type, ".txt"))
}
