# Filter coexpression matrices based on the proportion of cells in the dataset
# detectably expressing each gene. 
setwd("~/git/SCT-MoA")
library(methods)
library(dismay)
library(tidyverse)
library(magrittr)

# usage: filter-matrices.R <expr_dir> <coexpr_dir> <output_dir> optional: <idx>
args = commandArgs(trailingOnly = T)
if (length(args) < 3)
  stop("must provide expression, coexpression, and output directories")
expr_dir = args[1]
coexpr_dir = args[2]
output_dir = args[3]
if (!dir.exists(expr_dir))
  stop("expression directory does not exist: ", expr_dir)
if (!dir.exists(coexpr_dir))
  stop("coexpression directory does not exist: ", coexpr_dir)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)
# optionally, get file index (for multithreading)
fileIdx = NULL
if (length(args) >= 4)
  fileIdx = as.numeric(args[4])

# define thresholds
thresholds = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
# make sure output directories exist
threshold_dirs = file.path(output_dir, thresholds * 100)
for (threshold_dir in threshold_dirs)
  if (!dir.exists(threshold_dir))
    dir.create(threshold_dir, recursive = T)

# read files 
files = list.files(expr_dir, pattern = "*.gz", full.names = T)
if (!is.null(fileIdx))
  files = files[fileIdx]

# process each file in turn
for (file in files) {
  filename = gsub("\\.txt.*$", "", basename(file))
  message("filtering coexpression networks for file ", filename)
  dat = as.matrix(read.delim(file))
  
  # define gene lists
  gene_lists = map(thresholds, ~ names(which(colMeans(dat == 0) < .))) %>%
    setNames(thresholds)
  
  # filter each correlation matrix in turn
  coefs = dismay::metrics()
  for (coef in coefs) {
    message("  filtering ", coef, " matrix ...")
    
    # if all outputs exist, skip
    outputs = file.path(output_dir, thresholds * 100, paste0(
      filename, "-", coef, ".Rdata"))
    if (all(file.exists(outputs))) {
      message("skipping ", coef, ": all files exist")
      next
    }
    
    # load coexpression matrix
    coexpr_file = file.path(coexpr_dir, paste0(filename, "-", coef, ".Rdata"))
    load(coexpr_file, verbose = T) ## coexpr
    
    # filter based on proportion of dropouts
    for (threshold in thresholds) {
      genes = gene_lists[[as.character(threshold)]]
      filtered = coexpr[genes, genes]
      
      # write
      output_file = file.path(output_dir, threshold * 100, paste0(
        filename, "-", coef, ".Rdata"))
      save(filtered, file = output_file)
    }
  }
}

