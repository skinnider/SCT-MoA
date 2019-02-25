# Filter all GEO datasets.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/geo/preprocessing-functions.R")

# read human and mouse protein-coding genes
coding = read.delim("data/ensembl/protein_coding_genes.txt.gz")

# get all expression matrices
exprFiles = list.files("data/geo/processed", pattern = "*.gz", full.names = T)
for (exprFile in exprFiles) {
  accession = gsub("\\..*$", "", basename(exprFile))
  output_file = paste0("data/geo/filtered/", accession, ".txt")
  message("preprocessing file ", accession, " (", 
          which(exprFiles == exprFile), " of ", length(exprFiles), ") ...")
  
  # read data
  expr = as.matrix(read.delim(exprFile))
  
  # pre-filter to protein-coding genes only
  expr = expr[, colnames(expr) %in% coding$gene]
  
  # record rows and columns
  nrow1 = nrow(expr)
  ncol1 = ncol(expr)
  
  # replace NAs with zeroes
  expr[is.na(expr)] = 0
  
  # filter 1: remove samples that have no zeroes
  f1 = rowSums(expr == 0) > 0
  if (sum(!f1) > 0)
    message(".. removing ", sum(!f1), " \"cells\" with no zeroes!")
  expr = expr[f1, ]
  
  # filter 2: remove genes that are zero in >95% of samples
  f2 = colMeans(expr == 0) < 0.95
  expr = expr[, f2]
  
  # print filtering results
  message(".. writing filtered data with ", nrow(expr), " of ", nrow1, 
          " cells and ", ncol(expr), " of ", ncol1, " genes ...")
  
  # write filtered dataset
  write_and_gzip(expr, output_file)
}
