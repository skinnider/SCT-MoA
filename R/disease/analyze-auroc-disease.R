# Calculate accuracy of disease gene prediction in gene coexpression networks
# GSE75330_MFOL2 using a series of measures of correlation from
# single-cell RNA-seq datasets.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(flavin)
library(ontologyIndex)
library(EGAD)
library(dplyr)
library(data.table)
library(scales)
library(dismay)

# read arguments
# usage: analyze-auroc.R <directory> <outDirectory> <fileIdx>
args = commandArgs(trailingOnly = T)
if (length(args) < 3)
  stop("Must set input and output directories and input file index", call. = F)
directory = args[1]
outDirectory = args[2]
fileIdx = as.integer(args[3])
if (!dir.exists(directory))
  stop("input directory does not exist")
if (!dir.exists(outDirectory))
  dir.create(outDirectory, recursive = T)

# functions
source("R/Functions.R")

# check whether output file already exists
files = list.files(directory, full.names = T, pattern = "*\\.Rdata$")
## restrict disease analysis to:
## (1) Zeisel et al. 2018 mouse brain atlas
## (2) Vanlandewijck et al. 2018 vasculature atlas
files = files[grepl("-expr|GSE98816", files)]
file = files[fileIdx]
message(".. processing file ", file, " ...")
filename = basename(file)
coefs = dismay::metrics()
coef = last(coefs[sapply(coefs, function(c) grepl(c, filename))])
if (is.na(coef) | length(coef) == 0) {
  message(".. could not match file ", filename, " to coexpression method")
  next
}
idx = unlist(gregexpr(coef, filename))
dataset = substr(filename, 0, idx - 2)
output = file.path(outDirectory, paste0(dataset, "_", coef, ".txt.gz"))
if (file.exists(paste0(output, ".gz"))) {
  stop("output for file ", filename, " (dataset ", dataset,
       " with coefficient ", coef, ") already exists!")
} else {
  message(".. analyzing file ", filename, " (dataset ", dataset,
          " with coefficient ", coef, ") ...")
}

# load coexpression matrix
if (!file.exists(file)) {
  message("couldn't open file ", file)
  next
}
load(file, verbose = T) ## coexpr

# replace missing values with median (ZI kendall)
coexpr[is.na(coexpr)] = median(coexpr, na.rm = T)
# replace infinite values with the minimum (binomial)
coexpr[is.infinite(coexpr)] = min(coexpr, na.rm = T)

# scale metrics to [-1, 1] range
if (!dismay::is_bounded(coef)) {
  coexpr = scales::rescale(coexpr, to = c(-1, 1))
}

# read disease genes
dis = read.delim("data/disease/phenopedia/phenopedia.txt.gz")

# make EGAD coexpression network
## modified from EGAD::build_coexp_network to accommodate missing values
message("making EGAD coexpression network ...")
n = nrow(coexpr)
genes_in_mat = rownames(coexpr)
net = matrix(rank(coexpr, na.last = "keep", ties.method = "average"),
             nrow = n, ncol = n)
rownames(net) = colnames(net) = genes_in_mat
net = net / max(net, na.rm = T)
diag(net) = 1

# run EGAD
message("running guilt-by-association analysis ...")
disease_subset = dis %>%
  dplyr::filter(gene %in% genes_in_mat) %>%
  as.matrix()
disease_subset = disease_subset[, c(2, 1)]
if (nrow(disease_subset) == 0)
  stop("no disease terms annotated to network genes")
diseases = unique(disease_subset[, "disease"])
annotations = make_annotations(disease_subset, genes_in_mat, diseases)
gba = run_GBA(net, annotations, min = 1, max = ncol(net))

# record AUROCs
aurocs = gba[[1]][, "auc"]
# get number of proteins to which each term was annotated
all_terms = colSums(annotations)
n_terms = all_terms[names(aurocs)]
# write result
result = data.frame(dataset = dataset, coefficient = coef,
                    term = names(aurocs), auroc = aurocs,
                    n_proteins = n_terms,
                    pct_proteins = n_terms / length(genes_in_mat))
write.table(result, output, quote = F, sep = "\t", row.names = F)
system(paste("gzip --force", output))
