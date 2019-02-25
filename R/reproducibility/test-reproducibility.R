# For each pair of single-cell pancreas expression datasets of a given
# cell type, calculate the z score for the rank correlation between
# permuted coexpression matrices.
setwd("~/git/SCT-MoA")
library(methods)
library(magrittr)
library(dismay)
library(vegan)

# usage: test-reproducibility.R <input_dir> <index>
args = commandArgs(trailingOnly = T)
if (length(args) < 2)
  stop("must provide input directory and index (1-30)")
input_dir = args[1]
index = as.integer(args[2])
if (!dir.exists(input_dir))
  stop("input directory does not exist: ", input_dir)

# get combinations of datasets:
# GSE85241 - Muraro et al, Cell Syst 2016
# E-MTAB-5061 - Segersolpe et al, Cell Metab 2016
# GSE81547 - Enge et al, Cell 2017
# GSE81608 - Xin et al, Cell Metab 2016
# GSE84133 - Baron et al, Cell Syst 2016
datasets = c("GSE85241", "E-MTAB-5061", "GSE81547", "GSE81608",
             "GSE84133_human")
combinations = t(combn(datasets, 2))
patt = paste0(datasets, collapse = "|")
# combine with cell types
cell_types = c("alpha", "beta", "delta")
combinations = cbind(rbind(combinations, combinations, combinations),
                     rep(cell_types, each = nrow(combinations)))
# combine with coefficients 
coefs = dismay::metrics()
combinations = cbind(do.call(rbind, replicate(length(coefs), combinations, 
                                              simplify = F)) ,
                     rep(coefs, each = nrow(combinations)))

# get combination
combination = combinations[index, ]
cell_type = combination[3]
coef = combination[4]

# test whether file exists
output_dir = "data/reproducibility/pancreas"
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)
output = paste0(output_dir, "/", index, ".txt")
if (file.exists(output)) 
  stop("output file exists: ", output)
message("analyzing combination: ", paste(combination, collapse = " | "))

# seed RNG
set.seed(0)

# load matrices
files = file.path(input_dir, paste0(combinations[index, 1:2], "_", cell_type, 
                                    "-", coef, ".Rdata"))
load(files[1], verbose = F) ## filtered
coexpr1 = filtered
load(files[2], verbose = F) ## filtered
coexpr2 = filtered
rm(filtered)

# subset to common genes
common = intersect(rownames(coexpr1), rownames(coexpr2))
coexpr1 %<>% extract(common, common)
coexpr2 %<>% extract(common, common)

# replace missing values with the median
coexpr1 %<>% replace(!is.finite(.), median(., na.rm = T))
coexpr2 %<>% replace(!is.finite(.), median(., na.rm = T))

# run Mantel test
test = vegan::mantel(-coexpr1, -coexpr2, method = "spearman", 
                     permutations = 100)
stat = test$statistic
pval = test$signif
z_score = (test$statistic - mean(test$perm, na.rm = T)) / 
  sd(test$perm, na.rm = T)

# create data frame
df = data.frame(id1 = combination[1], id2 = combination[2],
                cell_type = combination[3], coefficient = combination[4],
                statistic = stat, z_score = z_score)

# write
write.table(df, output, quote = F, sep = "\t", row.names = F)
