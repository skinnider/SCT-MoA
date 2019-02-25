# Dataset: GSE63472
# Cells: single cells from mouse retinas
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")
library(data.table)

# read data
dat = as.data.frame(fread(
  "zcat < data/geo/raw/GSE63472_P14Retina_merged_digital_expression.txt.gz"))
# read cell type clusters
clust = read.delim("data/geo/soft/retina_clusteridentities.txt", header = F,
                   col.names = c("barcode", "cluster"))
# subset to valid cells
genes = dat[, 1]
dat = dat[, colnames(dat) %in% clust$barcode]

# analyze each of the largest clusters individually
clusterCounts = table(clust$cluster)
clusterCounts = clusterCounts[order(-clusterCounts)]
clusters = names(clusterCounts)[2:9]
for (cluster in clusters) {
  message(".. processing cluster ", cluster, " ...")
  barcodes = clust$barcode[clust$cluster == cluster]
  dat0 = dat[, colnames(dat) %in% barcodes]
  
  # map Entrez to ensembl
  genes0 = Hmisc::capitalize(tolower(genes))
  dat0 = map_genes(dat0, genes0, from = "SYMBOL", to = "ENSEMBL", 
                   db = org.Mm.eg.db)
  
  # convert to matrix and transpose
  genes1 = dat0[, 1]
  expr = t(dat0[, -1])
  colnames(expr) = genes1
  
  # print number of cells
  message("... writing ", nrow(expr), " cells in cluster ", cluster, " ...")
  
  # write
  output = paste0("data/geo/processed/GSE63472_c", cluster, ".txt")
  write_and_gzip(expr, output)
}
