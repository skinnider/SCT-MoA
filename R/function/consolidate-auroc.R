# Filter the complete set of AUROCs to (i) GO terms within the GO slim, and
# (ii) annotated to between 10 and 1,000 proteins in the network. 
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(ontologyIndex)
library(tidyverse)
library(magrittr)
source("R/Functions.R")

# usage: filter-auroc.R <dir>
args = commandArgs()
if (length(args) < 1)
  stop("must provide input directory")
dir = args[1]
if (!file.exists(dir))
  stop("input directory does not exist: ", dir)

# read GO slim
slim = get_ontology("data/go/goslim_generic.obo.gz")

# process files
all = list()
message("reading files in directory ", dir, " ...")
files = list.files(dir, pattern = "*.gz", full.names = T)
pb = progress::progress_bar$new(
  format = "file :what [:bar] :percent eta: :eta",
  clear = F, total = length(files), width = 100)
for (i in seq_len(length(files))) {
  file = files[i]
  dat = read.delim(file)
  filtered = dplyr::filter(dat, term %in% slim$id & 
                             dplyr::between(n_proteins, 10, 1000)) %>%
    dplyr::select(dataset, coefficient, term, auroc)
  all[[i]] = filtered
  # tick progress bar
  pb$tick(tokens = list(what = sprintf(
    paste0("%-", nchar(length(files)), "s"), i)))
}
combined = bind_rows(all)

# write 
write_and_gzip(combined, "data/function/auroc.txt")
