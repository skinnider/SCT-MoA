# Consolidate AUROCs for GO slim terms within networks at different levels of 
# non-zero gene expression filtering, and filter to terms
# annotated to between 10 and 1,000 proteins in the network. 
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(ontologyIndex)
library(tidyverse)
library(magrittr)
source("R/Functions.R")

# analyze each level of filtering separately 
dirs = list.dirs("data/function/auroc/filtered", recursive = F)
dat = data.frame()
for (dir in dirs) {
  threshold = basename(dir)
  message("consolidating networks at threshold: ", threshold, " ...")
  files = list.files(dir, pattern = "*.gz", full.names = T)
  df = map(files, read.delim) %>%
    bind_rows() %>%
    mutate(threshold = threshold) %>%
    dplyr::filter(dplyr::between(n_proteins, 10, 1000)) %>%
    dplyr::select(dataset, coefficient, threshold, term, auroc)
  dat %<>% rbind(df)
}

# write 
write_and_gzip(dat, "data/function/auroc_filtered.txt")
