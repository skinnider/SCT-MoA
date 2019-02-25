# Dataset: GSE98816
# Cells: single cells from the mouse brain vasculature
# Values: normalized counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim(
  "data/geo/raw/GSE98816/GSE98816_Brain_samples_normalized_counts_matrix.txt.gz")

# map symbols to Ensembl
symbols = dat[[1]]
map = AnnotationDbi::select(org.Mm.eg.db, keys = symbols, keytype = 'SYMBOL',
                            columns = 'ENSEMBL')
genes = map$ENSEMBL[match(symbols, map$SYMBOL)]

# convert to matrix and transpose
expr = t(dat[, -1])
colnames(expr) = genes
expr = expr[, !is.na(colnames(expr)) & !duplicated(colnames(expr))]

# read cell types
ct = read.csv("data/geo/raw/GSE98816/brain_3418cell_cluster_info.csv")
types = names(which(table(ct$cluster) >= 30))

# write each vascular cell type separately 
for (type in setdiff(types, c("AC", "OL"))) {
  cells = with(ct, cells[cluster == type])
  subset = expr[cells, ]
  output = paste0("data/geo/processed/GSE98816/GSE98816_", type, ".txt")
  write.and.gzip(subset, output)
  message("wrote cell type ", type, " with ", nrow(subset), " cells")
}
