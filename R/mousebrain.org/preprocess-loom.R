# Preprocess loom files.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)

# usage: preprocess-loom.R <input_dir> <output_dir>
args = commandArgs(trailingOnly = T)
if (length(args) < 2)
  stop("must provide input and output directories")
input_dir = args[1]
output_dir = args[2]
if (!dir.exists(input_dir))
  stop("input directory does not exist")
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)

# read protein-coding genes
coding = read.delim("data/ensembl/protein_coding_genes.txt.gz")

# process files one at a time
files = list.files(input_dir, pattern = "*.csv.gz", full.names = T)
for (file in files) {
  message("processing ", basename(file), " ...")
  
  # read data
  dat = read.csv(file, check.names = F)
  
  # filter empty genes
  genes = dat[[1]]
  expr = t(dat[, -1])
  colnames(expr) = genes
  expr = expr[rowSums(expr) > 0, colSums(expr) > 0]
  
  # filter to protein-coding genes
  expr = expr[, colnames(expr) %in% coding$gene]
  
  # filter to the top 5,000 genes
  present = colSums(expr > 0)
  ranks = rank(present, ties.method = 'random')
  keep = ranks > ncol(expr) - 5e3
  expr = expr[, keep]
  
  # write
  clean_name = chartr(' /', '__', gsub("\\.gz", "", basename(file)))
  output_file = file.path(output_dir, clean_name)
  write.table(expr, output_file, quote = F, sep = "\t", row.names = T)
  system(paste("gzip --force", output_file))
  message("  wrote file ", basename(file), " with ", nrow(expr), " cells and ",
          ncol(expr), " genes")
}
