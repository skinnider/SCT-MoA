# Consolidate information about scRNA-seq expression datasets into a single
# Excel file. 
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(openxlsx)

# read 1a
prop = read.csv("data/geo/datasets.csv")

# read 1b
drop = read.csv("data/geo/dropouts.csv", check.names = F)
colnames(drop)[1] = "Dataset"

# write to Excel file
write.xlsx(list(
  "Supplementary Table 1a" = prop, "Supplementary Table 1b" = drop), 
  "data/tables/Table_S1.xlsx")
