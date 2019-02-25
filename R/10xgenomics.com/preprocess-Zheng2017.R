# Preprocess ten human bead-enriched PBMC subpopulations described by
# Zheng et al., Nature Communications 2017, and obtained from
# support.10xgenomics.com.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# usage: preprocess-Zheng2017.R <input_dir> <output_dir>
args = commandArgs(trailingOnly = T)
if (length(args) < 2)
  stop("must provide input and output directories")
input_dir = args[1]
output_dir = args[2]
if (!dir.exists(input_dir))
  stop("input directory does not exist")
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)

# get list of files from Zhang et al., 2017. downloaded from:
# https://support.10xgenomics.com/single-cell-gene-expression/datasets
files = list.files(input_dir, pattern = '.tar.gz$', full.names = T)
for (file in files) {
  message(".. working on file ", which(files == file), " of ", length(files), 
          ": ", basename(file))
  
  # get file information
  cell_type = gsub('_filtered_gene_bc_matrices.tar.gz', '', basename(file))
  
  # untar
  matrix_dir = file.path(input_dir, 'mat', cell_type)
  if (!dir.exists(matrix_dir))
    dir.create(matrix_dir, recursive = T)
  untar(file, exdir = path.expand(matrix_dir))
  
  # create data matrix
  mat_files = list.files(matrix_dir, recursive = T, full.names = T)
  mat = as.matrix(Matrix::readMM(mat_files[grep(".mtx", mat_files)]))
  genes = read.delim(mat_files[grep("genes", mat_files)], header = F)[, 1]
  mat = t(mat)
  colnames(mat) = genes

  # write
  write_and_gzip(mat, file.path(
    output_dir, paste0('Zheng_', cell_type, '.txt')))
}
