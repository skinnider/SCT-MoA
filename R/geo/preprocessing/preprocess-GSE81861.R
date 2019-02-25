# Dataset: GSE81861
# Cells: single cells from seven cell lines
# Values: FPKM
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# load data 
dat = read.csv("data/geo/raw/GSE81861/GSE81861_Cell_Line_FPKM.csv.gz")

# fix gene labels
genes = gsub("^.*_|\\..*$", "", dat[[1]])

# convert to matrix
expr = t(as.matrix(dat[, -1]))
colnames(expr) = genes

# fix cell lines
split = strsplit(rownames(expr), "__")
lines = map_chr(split, 2)
rownames(expr) = lines

# write 
filepath = "data/geo/processed/GSE81861.txt"
write.table(expr, filepath, sep = "\t", quote = F, row.names = T)
system(paste("gzip --force", filepath))

# also filter here, as this dataset is not used in coexpression analysis
# read human and mouse protein-coding genes
coding = read.delim("data/ensembl/protein_coding_genes.txt.gz")
# protein-coding genes only
expr = expr[, colnames(expr) %in% coding$gene]
# remove genes that are zero in 80% of samples
f2 = colMeans(expr == 0) < 0.8
expr = expr[, f2]
# write 
filepath = "data/geo/filtered/GSE81861.txt"
write.table(expr, filepath, sep = "\t", quote = F, row.names = T)
system(paste("gzip --force", filepath))
