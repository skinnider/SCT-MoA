# Consolidate reproducibility test results into a single file. 
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# get files
files = list.files("data/reproducibility/pancreas", pattern = "*.txt", 
                   full.names = T)
# consolidate
dat = map(files, read.delim) %>%
  bind_rows()
# write
write.table(dat, "data/reproducibility/pancreas.txt", quote = F,
            sep = "\t", row.names = F)
