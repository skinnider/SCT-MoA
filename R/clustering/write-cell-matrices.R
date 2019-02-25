# Construct cell correlation matrices for RCA cell lines dataset.
setwd("~/git/SCT-MoA")
library(methods)
library(dismay)

# read file
file = "data/geo/filtered/GSE81861.txt.gz"
dat = read.delim(file, row.names = NULL)
names = dat[[1]]
cell = t(as.matrix(dat[, -1]))
colnames(cell) = paste0(names, '-', seq_len(ncol(cell)))

# calculate each correlation
coefs = dismay::metrics()
for (coef in coefs) {
  message("  calculating ", coef, " matrix ...")
  # create correlation matrix 
  output = file.path("data/clustering/matrices", paste0(
    "GSE81861-", coef, ".Rdata"))
  if (!file.exists(output)) {
    coexpr = dismay::dismay(cell, metric = coef)
    save(coexpr, file = output)
  }
}
