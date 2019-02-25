# Dataset: GSE60768
# Cells: mouse embryonic stem cells and neural stem cells 
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.csv("data/geo/raw/GSE60768_TotalGeneCounts.csv.gz")
genes = dat[, 1]

# process NSCs and mESCs separately
patterns = c("NSCs" = "_NSC", "mESCs" = "_2i|_ESC")
for (cellType in names(patterns)) {
  patt = patterns[cellType]
  
  # subset to relevant cells
  dat0 = dat[, grepl(patt, colnames(dat))]
  
  # convert to matrix and transpose
  expr = t(dat0)
  colnames(expr) = genes
  
  # write
  write_and_gzip(expr, paste0("data/geo/processed/GSE60768_", cellType, ".txt"))
}
