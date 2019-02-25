# Consolidate all network overlap data into a single file and gzip it. 
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)

# list files
files = list.files("data/networks/overlap", pattern = "*.txt", full.names = T)
# read data
dat = map(files, read.delim) %>%
  bind_rows()
# write
output = "data/networks/overlap.txt"
write.table(dat, output, quote = F, sep = "\t", row.names = F)
system(paste("gzip --force", output))
