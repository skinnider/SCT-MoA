# Dataset: GSE87544
# Cells: single cells from the mouse hypothalamus
# Values: TPM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# load data 
load(
  "data/geo/raw/GSE87544/GSE87544_1443737Cells.Expresssion.Matrix.log_tpm+1_.renamed.RData", 
  verbose = T)
## rename
dat = Expresssion_Matrix_unfiltered
rm(Expresssion_Matrix_unfiltered)

# map symbols to Ensembl
symbols = rownames(dat)
map = AnnotationDbi::select(org.Mm.eg.db, keys = symbols, keytype = 'SYMBOL',
                            columns = 'ENSEMBL')
genes = map$ENSEMBL[match(symbols, map$SYMBOL)]

# convert to matrix
expr = t(dat)
colnames(expr) = genes
expr = expr[, !is.na(colnames(expr)) & !duplicated(colnames(expr))]

# read cell types
types = read.csv(
  "data/geo/raw/GSE87544/GSE87544_1443737Cells.SVM.cluster.identity.renamed.csv.gz")
# write each separately
major_types = names(which(table(types$SVM_clusterID) >= 100))
major_types %<>% setdiff("zothers")
for (type in major_types) {
  cells = types$Cell_ID_temp[types$SVM_clusterID == type]
  subset = expr[cells, ]
  output = paste0("data/geo/processed/GSE87544_", type, ".txt")
  write_and_gzip(subset, output)
  message("wrote cell type ", type, " with ", nrow(subset), " cells")
}
