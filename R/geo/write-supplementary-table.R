# Write information about the number of cells, number of genes, and proportion
# of dropouts in each dataset.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# first, collect all files
files = list.files("data/geo/filtered", pattern = "gz") %>% 
  gsub("\\.txt.*$", "", .)
loom = list.files("data/loom/raw", pattern = "gz") %>% 
  gsub("\\.csv.*$", "", .)
files %<>% c(loom)
df = data.frame(dataset = gsub("-expression", "", files),
                cells = NA, genes = NA, dropouts = NA)

# for GEO datasets, calculate # cells, # genes, % dropouts
for (file_idx in seq_along(files)) {
  file = files[file_idx]
  message("processing file ", file_idx, " of ", length(files), ": ", file)
  mat = as.matrix(read.delim(paste0("data/geo/filtered/", file, ".txt.gz")))
  keep = colMeans(mat == 0) < 0.8
  mat %<>% extract(, keep)
  cells = nrow(mat)
  genes = ncol(mat)
  dropouts = mean(mat == 0, na.rm = T)
  df[file_idx, 2:4] = c(cells, genes, dropouts)
}

# get the same data for loom matrices
loom_files = list.files("data/loom/raw", pattern = "*.gz", full.names = T)
loom = map(loom_files, ~ as.matrix(read.delim(.)))
cells = map_int(loom, nrow)
dropouts = map_dbl(loom, ~ mean(. == 0))
df[match(gsub("-expr.*$", "", basename(loom_files)), df$dataset), 2:4] = 
  data.frame(cells = cells, genes = NA, dropouts = dropouts)

# write to CSV
write.csv(df, "data/geo/datasets.csv", quote = F, row.names = F)

## manually merge in curated data on PMID, species, and protocol
