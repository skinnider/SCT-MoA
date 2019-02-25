# Dataset: GSE99058
# Cells: single mouse astrocytes
# Values: normalized counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim(
  "data/geo/raw/GSE99058_Brain_AC_250samples_normalized_counts_matrix.txt.gz")

# map symbols to Ensembl
symbols = dat[[1]]
map = AnnotationDbi::select(org.Mm.eg.db, keys = symbols, keytype = 'SYMBOL',
                            columns = 'ENSEMBL')
genes = map$ENSEMBL[match(symbols, map$SYMBOL)]

# convert to matrix and transpose
expr = t(dat[, -1])
colnames(expr) = genes
expr = expr[, !is.na(colnames(expr)) & !duplicated(colnames(expr))]

# write
output = paste0("data/geo/processed/GSE99058.txt")
write_and_gzip(expr, output)
