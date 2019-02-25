# Dataset: GSE81547
# Cells: single cells from human pancreatic islets
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")
library(tidyverse)

# read expression data 
expr_files = list.files("data/geo/raw/GSE81547", pattern = "*.csv.gz", 
                        full.names = T)
expr = map(expr_files, ~ read.delim(., header = F)) %>%
  bind_cols()

# map symbols to Ensembl
symbols = expr[[1]]
map = AnnotationDbi::select(org.Hs.eg.db, keys = symbols, keytype = 'SYMBOL',
                            columns = 'ENSEMBL')
genes = map$ENSEMBL[match(symbols, map$SYMBOL)]

# convert to count matrix and map to cell names
expr_idxs = map_lgl(expr, is.integer)
expr_accns = gsub("_.*$", "", basename(expr_files))
expr = t(expr[, expr_idxs])
rownames(expr) = expr_accns
colnames(expr) = genes
expr = expr[, !is.na(colnames(expr)) & !duplicated(colnames(expr))]

# read cell types from SOFT file
lines = readLines("data/geo/raw/GSE81547/GSE81547_family.soft.gz")
accns = lines[grepl("Sample_geo_accession", lines)]
types = lines[grepl("inferred_cell_type", lines)]

# map cell types to accessions
types = gsub("^.*\\: ", "", types) %>%
  setNames(gsub("^.*= ", "", accns))

# split up by cell types
cell_types = setdiff(unique(types), 'unsure')
for (type in cell_types) {
  cells = names(which(types == type))
  subset = expr[rownames(expr) %in% cells,]
  output_file = paste0("data/geo/processed/GSE81547_", type, ".txt")
  write_and_gzip(subset, output_file)
  message("wrote cell type `", type, "` with ", nrow(subset), " cells")
}
