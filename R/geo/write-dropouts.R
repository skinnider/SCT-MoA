# Write the total proportion of dropouts in each dataset at various gene
# filtering thresholds.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# process GEO, 10X, and loom separately
files_geo = list.files("data/geo/processed", pattern = "*.gz$", full.names = T)
files_10x = list.files("data/10xgenomics.com/processed", pattern = "*.gz$", 
                       full.names = T)
files_loom = list.files("data/loom/filtered", pattern = "*.gz$", full.names = T)

# define thresholds
thresholds = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, NA)
threshold_names = paste0("dropouts, ", 100 * thresholds, "%")
threshold_names[length(thresholds)] = "dropouts, original dataset"

# process GEO files
df = data.frame()
for (file_idx in seq_along(files_geo)) {
  file = files_geo[file_idx]
  if (grepl("GSE81861", file))
    next
  filename = gsub("\\..*$", "", basename(file))
  message("processing GEO file ", file_idx, " of ", length(files_geo), ": ", 
          filename)
  expr = as.matrix(read.delim(file))
  dropouts = map_dbl(thresholds, ~ ifelse(
    is.na(.), 
    mean(expr == 0),
    mean(expr[, colMeans(expr == 0) < .] == 0)))
  row = data.frame(dataset = filename)
  for (idx in seq_along(threshold_names))
    row[[threshold_names[idx]]] = dropouts[idx]
  df %<>% rbind(row)
}

# process 10X files
for (file_idx in seq_along(files_10x)) {
  file = files_10x[file_idx]
  filename = gsub("\\..*$", "", basename(file))
  message("processing 10X file ", file_idx, " of ", length(files_10x), ": ", 
          filename)
  expr = as.matrix(read.delim(file))
  dropouts = map_dbl(thresholds, ~ ifelse(
    is.na(.), 
    mean(expr == 0),
    mean(expr[, colMeans(expr == 0) < .] == 0)))
  row = data.frame(dataset = filename)
  for (idx in seq_along(threshold_names))
    row[[threshold_names[idx]]] = dropouts[idx]
  df %<>% rbind(row)
}

# process loom files
for (file_idx in seq_along(files_loom)) {
  file = files_loom[file_idx]
  filename = gsub("-expr.*$", "", basename(file))
  message("processing loom file ", file_idx, " of ", length(files_loom), ": ", 
          filename)
  expr = as.matrix(read.delim(file))
  dropouts = mean(expr == 0)
  row = data.frame(dataset = filename)
  for (idx in seq_along(threshold_names))
    row[[threshold_names[idx]]] = NA
  row[["dropouts, original dataset"]] = dropouts
  df %<>% rbind(row)
}

# write
write.csv(df, "data/geo/dropouts.csv", row.names = F)
