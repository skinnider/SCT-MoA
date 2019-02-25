# Dataset: E-MTAB-2805
# Cells: mouse ESCs staged for the G1, G2M, and S cell-cycle phases
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")

# combine cells at each phase (low n)
stages = c("G1", "G2M", "S")
stageFiles = paste0("data/geo/raw/E-MTAB-2805/", stages, 
                    "_singlecells_counts.txt.gz")
dats = lapply(stageFiles, read.delim)
dat = do.call(data.frame, dats)

# extract genes
genes = dat[, 1]

# remove genes, transcripts, symbols, and lengths
dat = dat[, grepl("count", colnames(dat))]

# convert to matrix and transpose
expr = t(dat)
colnames(expr) = genes

# write expression
write_and_gzip(expr, "data/geo/processed/E-MTAB-2805.txt")
