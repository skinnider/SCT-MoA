# Dataset: GSE85241
# Cells: single cells from human pancreatic islets
# Values: normalized counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim(
  "data/geo/raw/GSE85241/GSE85241_cellsystems_dataset_4donors_updated.txt.gz")

# map symbols to Ensembl
symbols = gsub("__.*$", "", rownames(dat))
map = AnnotationDbi::select(org.Hs.eg.db, keys = symbols, keytype = 'SYMBOL',
                            columns = 'ENSEMBL')
genes = map$ENSEMBL[match(symbols, map$SYMBOL)]

# convert to matrix and transpose
expr = t(dat)
colnames(expr) = genes
expr = expr[, !is.na(colnames(expr)) & !duplicated(colnames(expr))]

# split up by cell types
types = read.delim("data/geo/raw/GSE85241/cell_type_annotation_Cels2016.txt")
cell_types = setdiff(unique(types$x), c("epsilon", "unclear"))
for (type in cell_types) {
  cells = rownames(types)[types$x == type]
  subset = expr[rownames(expr) %in% cells,]
  output_file = paste0("data/geo/processed/GSE85241_", type, ".txt")
  write_and_gzip(subset, output_file)
  message("wrote cell type `", type, "` with ", nrow(subset), " cells")
}
