# Consolidate functional coherence and network overlap analyses into 
# supplementary datasets.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(openxlsx)
source("R/functions.R")

# read AUC
auc = read.delim("data/function/auroc.txt.gz") %>%
  mutate(dataset = gsub("-expr.*$", "", dataset),
         coefficient = clean_metric(coefficient)) %>%
  set_colnames(c("Dataset", "Measure of association", "GO term", "AUROC"))

# read network overlap
ovr = read.delim("data/networks/overlap.txt.gz") %>%
  mutate(dataset = gsub("-expr.*$", "", dataset),
         coefficient = clean_metric(coefficient)) %>% 
  dplyr::select(dataset, coefficient, network, cutoff, obs, 
                rnd_mean, rnd_median, rnd_sd) %>%
  set_colnames(c("Dataset", "Measure of association", "Network", 
                 "Number of edges", "Observed overlap", 
                 "Randomized overlap, mean", "Randomized overlap, median",
                 "Randomized overlap, s.d."))

# write to Excel files
write.xlsx(list("Supplementary Data 1" = auc), "data/tables/Data_S1.xlsx")
write.xlsx(list("Supplementary Data 2" = ovr), "data/tables/Data_S2.xlsx")
