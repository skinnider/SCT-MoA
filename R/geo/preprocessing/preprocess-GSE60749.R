# Dataset: GSE60749
# Cells: mouse embryonic stem cells and neural precursor cells
# Values: TPM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# process mESCs and NPCs separately
base = "data/geo/raw/GSE60749/GSE60749_RnaSeq_single_cell_"
files = c("mESCs" = paste0(base, "mESC", "_TPM.csv.gz"),
          "NPCs" = paste0(base, "NPC", "_TPM.csv.gz"))
for (cellType in names(files)) {
  file = files[cellType]
  
  # read data
  dat = read.csv(file)
  genes = dat[, 1]
  dat = dat[, -1]
  
  # map to ensembl
  dat = map_genes(dat, genes, from = "SYMBOL", to = "ENSEMBL", 
                  db = org.Mm.eg.db)
  
  # convert to matrix and transpose
  genes = dat[, 1]
  expr = t(dat[, -1])
  colnames(expr) = genes
  
  # write
  write_and_gzip(expr, paste0("data/geo/processed/GSE60749_", cellType, ".txt"))
}
