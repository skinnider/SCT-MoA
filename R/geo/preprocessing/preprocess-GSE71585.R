# Dataset: GSE71585
# Cells: single cells from the primary visual cortex
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.csv("data/geo/raw/GSE71585/GSE71585_RefSeq_counts.csv.gz")

# map symbols to Ensembl
symbols = dat[[1]]
map = AnnotationDbi::select(org.Mm.eg.db, keys = symbols, keytype = 'SYMBOL',
                            columns = 'ENSEMBL')
genes = map$ENSEMBL[match(symbols, map$SYMBOL)]

# convert to matrix and transpose
expr = t(dat[, -1])
colnames(expr) = genes
expr = expr[, !is.na(colnames(expr)) & !duplicated(colnames(expr))]

# read clustering results
clust = read.csv("data/geo/raw/GSE71585/GSE71585_Clustering_Results.csv.gz")
# subset to cells in matrix
clust = clust[clust$sample_title %in% rownames(expr), ]
# remove intermediate or unclassified cells
clust = clust[clust$core_intermediate == 'core', ]

# write finer cell types
minor_types = names(which(table(clust$primary_type) >= 30))
for (type in minor_types) {
  cells = clust$sample_title[clust$primary_type == type]
  subset = expr[cells,]
  message(".. writing cell type ", type, " with ", nrow(subset), " cells ...")
  
  # write
  cleanType = chartr(' /', '--', gsub("-", "", type))
  output = paste0("data/geo/processed/GSE71585_l2_", cleanType, ".txt")
  write_and_gzip(subset, output)
}
