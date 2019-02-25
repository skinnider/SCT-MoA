# Construct gene correlation matrices for a given dataset, using 17
# measures of association.
setwd("~/git/SCT-MoA")
library(methods)
library(dismay)

# usage: write-matrices.R <directory> <outDirectory> optional: <idx>
args = commandArgs(trailingOnly = T)
if (length(args) == 0)
  stop("Must set input file directory", call. = F)
directory = args[1]
if (length(args) == 1)
  stop("Must set output file directory", call. = F)
outDirectory = args[2]
# optionally, get file index (for multithreading)
fileIdx = NULL
if (length(args) >= 3) 
  fileIdx = as.numeric(args[3])

# read files 
files = list.files(directory, pattern = "*.gz", full.names = T)
if (!is.null(fileIdx))
  files = files[fileIdx]

# process each file in turn
for (file in files) {
  filename = gsub("\\.txt.*$", "", basename(file))
  message("creating coexpression networks for file ", filename)
  dat = as.matrix(read.delim(file))

  # calculate each correlation
  coefs = dismay::metrics()
  for (coef in coefs) {
    message("  calculating ", coef, " matrix ...")
    # create correlation matrix 
    output = file.path(outDirectory, paste0(filename, "-", coef, ".Rdata"))
    if (!file.exists(output)) {
      coexpr = dismay::dismay(dat, metric = coef)
      save(coexpr, file = output)
    }
  }
}
