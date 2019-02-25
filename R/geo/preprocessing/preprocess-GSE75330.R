# Dataset: GSE75330
# Cells: single oligodendrocytes from mouse brain regions
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim(
  "data/geo/raw/GSE75330/GSE75330_Marques_et_al_mol_counts2.tab.gz", 
  check.names = F)

# read SOFT files
lines = readLines("data/geo/raw/GSE75330/GSE75330_family.soft.gz")
samples = lines[grepl("Sample_title", lines)]
types = lines[grepl("inferred cell type", lines)]
cell_types = setNames(gsub("^.*\\: ", "", types), gsub("^.*= ", "", samples))

# match symbols to Ensembl
symbols = dat$cellid
map = AnnotationDbi::select(org.Mm.eg.db, keys = symbols, keytype = 'SYMBOL',
                            columns = 'ENSEMBL')
genes = map$ENSEMBL[match(symbols, map$SYMBOL)]

# convert to matrix
expr = t(dat[, -1])
colnames(expr) = genes
expr = expr[, !is.na(colnames(expr)) & !duplicated(colnames(expr))]

# split up by cell types
for (type in unique(cell_types)) {
  cells = names(which(cell_types == type))
  subset = expr[rownames(expr) %in% cells,]
  output_file = paste0("data/geo/processed/GSE75330_", type, ".txt")
  write_and_gzip(subset, output_file)
  message("wrote cell type `", type, "` with ", nrow(subset), " cells")
}
