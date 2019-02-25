# Dataset: GSE84133
# Cells: single cells from human and mouse pancreatic islets
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")
library(tidyverse)

# read orthologs
orthologs = read.delim("data/orthologs/human-mouse-orthologs.txt", 
                       col.names = c("human", "mouse"))
# subset to one-to-one
oto = orthologs[with(orthologs, !duplicated(human) & !duplicated(mouse) & 
                       human != "" & mouse != ""), ]

# list human and mouse files
expr_files = list.files("data/geo/raw/GSE84133", pattern = "*.csv.gz",
                        full.names = T)
hsa_files = expr_files[grepl("human", expr_files)]
mmu_files = expr_files[grepl("mouse", expr_files)]

# read expression
hsa = map(hsa_files, ~ read.csv(., check.names = F)) %>% 
  bind_rows()
mmu = map(mmu_files, ~ read.csv(., check.names = F)) %>% 
  bind_rows()

# map symbols to Ensembl
hsa_symbols = colnames(hsa)[map_lgl(hsa, is.numeric)]
mmu_symbols = colnames(mmu)[map_lgl(mmu, is.numeric)]
hsa_map = AnnotationDbi::select(org.Hs.eg.db, keys = hsa_symbols, 
                                keytype = 'SYMBOL', columns = 'ENSEMBL')
mmu_map = AnnotationDbi::select(org.Mm.eg.db, keys = mmu_symbols, 
                                keytype = 'ALIAS', columns = 'ENSEMBL')
hsa_genes = hsa_map$ENSEMBL[match(hsa_symbols, hsa_map$SYMBOL)]
mmu_genes = mmu_map$ENSEMBL[match(mmu_symbols, mmu_map$ALIAS)]
# match mouse to orthologs
hsa_orthologs = oto$human[match(mmu_genes, oto$mouse)]

# split up by cell types
cell_types = unique(c(hsa$assigned_cluster, mmu$assigned_cluster))
exclude = c("t_cell", "T_cell", "immune_other", "macrophage", "epsilon", "mast",
            "schwann", "activated_stellate", "quiescent_stellate", "B_cell")
cell_types %<>% setdiff(exclude)
for (type in cell_types) {
  # write human 
  hsa_subset = hsa[hsa$assigned_cluster == type, ]
  hsa_expr = as.matrix(hsa_subset[, map_lgl(hsa_subset, is.numeric)])
  colnames(hsa_expr) = hsa_genes
  hsa_expr = hsa_expr[, !is.na(colnames(hsa_expr)) & 
                        !duplicated(colnames(hsa_expr))]
  output_file = paste0("data/geo/processed/GSE84133_human_", type, ".txt")
  write_and_gzip(hsa_expr, output_file)
  message("wrote cell type `", type, "` with ", nrow(hsa_expr),
          " human cells")
  
  # write mouse 
  mmu_subset = mmu[mmu$assigned_cluster == type, ]
  if (nrow(mmu_subset) == 0)
    next
  mmu_expr = as.matrix(mmu_subset[, map_lgl(mmu_subset, is.numeric)])
  colnames(mmu_expr) = hsa_orthologs
  mmu_expr = mmu_expr[, !is.na(colnames(mmu_expr)) & 
                        !duplicated(colnames(mmu_expr))]
  output_file = paste0("data/geo/processed/GSE84133_mouse_", type, ".txt")
  write_and_gzip(mmu_expr, output_file)
  message("wrote cell type `", type, "` with ", nrow(mmu_expr),
          " mouse cells")
}
