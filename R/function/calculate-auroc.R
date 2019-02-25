# Analyze the functional connectivity of gene coexpression networks constructed
# using a series of measures of correlation for single-cell RNA-seq datasets.
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
# usage: calculate-auroc.R <directory> <outDirectory> <idx>
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
source("R/functions.R")

# check whether output file already exists
files = list.files(directory, full.names = T, pattern = "*\\.Rdata$")
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
if (file.exists(output)) {
  stop("output for file ", filename, " (dataset ", dataset, 
       " with coefficient ", coef, ") already exists!")
} else {
  message(".. analyzing file ", filename, " (dataset ", dataset, 
          " with coefficient ", coef, ") ...")
}

# read protein-coding genes
coding = read.delim("data/ensembl/protein_coding_genes.txt.gz")
get_species = function(vector) {
  subset = coding %>% dplyr::filter(gene %in% vector)
  names(which.max(table(subset$source)))
}

# load coexpression matrix
if (!file.exists(file)) {
  message("couldn't open file ", file) 
  next      
}
load(file, verbose = T) ## coexpr

# replace missing values with median (ZI kendall)
if (any(is.na(coexpr)))
  coexpr[is.na(coexpr)] = median(coexpr, na.rm = T)
# replace infinite values with the minimum (binomial)
if (any(is.infinite(coexpr)))
  coexpr[is.infinite(coexpr)] = min(coexpr, na.rm = T)
  
# scale metrics to [-1, 1] range
if (!dismay::is_bounded(coef)) {
  coexpr = scales::rescale(coexpr, to = c(-1, 1))
}

# read GO terms
ontology = get_ontology("data/go/go-basic.obo.gz")
human = read_gpa("data/go/goa_human.gpa.gz", accession = "ENSEMBL", 
                 database = org.Hs.eg.db, ontology = ontology, propagate = T)
mouse = read_gpa("data/go/goa_mouse.gpa.gz", accession = "ENSEMBL", 
                 database = org.Mm.eg.db, ontology = ontology, propagate = T)
human = filter_breadth(human, min = 1, max = 1e4)
mouse = filter_breadth(mouse, min = 1, max = 1e4)

# detect species
species = get_species(colnames(coexpr))
# load species-specific GO file
go = get(species)

# make EGAD coexpression network
## modified from EGAD::build_coexp_network to accommodate missing values
message("making EGAD coexpression network ...")
n = nrow(coexpr)
genes = rownames(coexpr)
net = matrix(rank(coexpr, na.last = "keep", ties.method = "average"), 
             nrow = n, ncol = n)
rownames(net) = colnames(net) = genes
net = net / max(net, na.rm = T)
diag(net) = 1

# run EGAD
message("running guilt-by-association analysis ...")
go_subset = go %>%
  dplyr::filter(ENSEMBL %in% genes) %>%
  dplyr::select(ENSEMBL, GO.ID) %>%
  as.matrix()
if (nrow(go_subset) == 0)
  stop("no GO terms annotated to network genes. check right GO file")
go_terms = unique(go_subset[, "GO.ID"])
annotations = make_annotations(go_subset, genes, go_terms)
gba = run_GBA(net, annotations, min = 0, max = 1e4)

# record AUROCs 
aurocs = gba[[1]][, "auc"]
# get number of proteins to which each term was annotated
all_terms = colSums(annotations)
n_terms = all_terms[names(aurocs)]
# write result
result = data.frame(dataset = dataset, coefficient = coef, 
                    term = names(aurocs), auroc = aurocs, 
                    n_proteins = n_terms, 
                    pct_proteins = n_terms / length(genes))
write.table(result, output, quote = F, sep = "\t", row.names = F)
system(paste("gzip --force", output))
