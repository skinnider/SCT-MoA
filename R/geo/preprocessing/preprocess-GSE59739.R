# Dataset: GSE59739
# Cells: single sensory neurons from the dorsal root ganglion
# Values: RPM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.csv("data/geo/raw/GSE59739_Table1.csv.gz")

# subset to cells
symbols = dat[[1]]
expr = dat[, which(dat[9,] == 'cell')]
# subset to neuronal cells
expr = expr[, which(expr[6,] %in% c('NF', 'NP', 'PEP', 'TH'))]
# remove metadata
expr = expr[11:nrow(expr), ]
expr = as.matrix(expr)
mode(expr) = "numeric"

# map symbols to Ensembl
map = AnnotationDbi::select(org.Mm.eg.db, keys = symbols[-c(1:10)], 
                            keytype = 'SYMBOL', columns = 'ENSEMBL')
genes = map$ENSEMBL[match(symbols[-c(1:10)], map$SYMBOL)]
rownames(expr) = genes

# transpose
expr = t(expr)
expr = expr[, !is.na(colnames(expr)) & !duplicated(colnames(expr))]

# write
write_and_gzip(expr, "data/geo/processed/GSE59739.txt")
