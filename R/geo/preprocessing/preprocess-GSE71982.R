# Dataset: GSE71982
# Cells: single cells from the utricular and cochlear sensory epithelia
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# read data
dat = read.delim("data/geo/raw/GSE71982/GSE71982_RSEM_Counts_Matrix.txt.gz")

# remove gene names from the data frame 
genes = dat[, 1]
dat = dat[, -1]

# fix genes
genes = gsub("\\\"", "", genes)

# map to ensembl
dat = map_genes(dat, genes, from = "SYMBOL", to = "ENSEMBL", db = org.Mm.eg.db)

# convert to matrix and transpose
genes = dat[, 1]
expr = t(dat[, -1])
colnames(expr) = genes

# subset to analyzed cells
pheno = read.delim(
  "data/geo/raw/GSE71982/GSE71982_P1_Utricle_PhenoData.txt.gz")
expr = expr[pheno$Short_Name, ]

# write
write_and_gzip(expr, "data/geo/processed/GSE71982.txt")
