# Dataset: E-MTAB-5061
# Cells: single cells from human pancreatic islets
# Values: RPKM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim(
  "data/geo/raw/E-MTAB-5061/pancreas_refseq_rpkm_3514sc.txt.gz")

# map symbols to Ensembl
symbols = dat[, 1]
map = AnnotationDbi::select(org.Hs.eg.db, keys = symbols, keytype = 'SYMBOL',
                            columns = 'ENSEMBL')
genes = map$ENSEMBL[match(symbols, map$SYMBOL)]

# convert to matrix and transpose
expr = t(dat[, -c(1:2)])
colnames(expr) = genes
expr = expr[, !is.na(colnames(expr)) & !duplicated(colnames(expr))]

# split up by cell types
sdrf = read.delim("data/geo/raw/E-MTAB-5061/E-MTAB-5061.sdrf.txt.gz")
cell_types = setdiff(unique(sdrf$Factor.Value.cell.type.), c(
  "co-expression cell", "epsilon cell", "mast cell", "MHC class II cell",
  "not applicable", "unclassified cell", "unclassified endocrine cell",
  "endothelial cell"))
for (type in cell_types) {
  cells = sdrf$Source.Name[sdrf$Factor.Value.cell.type. == type]
  subset = expr[rownames(expr) %in% cells,]
  output_file = paste0("data/geo/processed/E-MTAB-5061_",
                       gsub(" .*$", "", type), ".txt")
  write_and_gzip(subset, output_file)
  message("wrote cell type `", type, "` with ", nrow(subset), " cells")
}
