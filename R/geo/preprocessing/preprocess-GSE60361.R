# Dataset: GSE60361
# Cells: single cells from two regions of the mouse cerebral cortex
# Values: count
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim("data/geo/raw/GSE60361/GSE60361_C1-3005-Expression.txt.gz")

# remove gene names from the data frame 
symbols = dat[, 1]
dat = dat[, -1]

# read sample table
samples = read.delim("data/geo/raw/GSE60361/GSE60361.txt")

# write minor cell types, regardless of region
minor_types = setdiff(names(which(table(samples$level2class) >= 30)), "(none)")
# map symbol to ensembl
map = AnnotationDbi::select(org.Mm.eg.db, keys = symbols, keytype = 'SYMBOL',
                            columns = 'ENSEMBL')
genes = map$ENSEMBL[match(symbols, map$SYMBOL)]
# convert to matrix
expr = t(dat)
colnames(expr) = genes
expr = expr[, !is.na(colnames(expr)) & !duplicated(colnames(expr))]
# split up by cell types
for (type in minor_types) {
  cells = samples$cell_id[samples$level2class == type]
  subset = expr[rownames(expr) %in% paste0("X", cells),]
  output_file = paste0("data/geo/processed/GSE60631_l2_", type, ".txt")
  write_and_gzip(subset, output_file)
  message("wrote cell type `", type, "` with ", nrow(subset), " cells")
}
