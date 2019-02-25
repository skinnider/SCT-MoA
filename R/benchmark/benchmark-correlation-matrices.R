# Benchmark the wall clock time it takes to construct gene coexpression
# matrices with each measure of association. 
setwd("~/git/SCT-MoA")
library(methods)
library(dismay)
library(microbenchmark)
library(magrittr)

# usage: benchmark-correlation-matrices.R <directory> <fileIdx>
args = commandArgs(trailingOnly = T)
if (length(args) == 0)
  stop("Must set input file directory", call. = F)
directory = args[1]
if (length(args) >= 2)
  fileIdx = as.integer(args[2])

# read files 
files = list.files(directory, pattern = "*.gz", full.names = T)
if (!is.null(fileIdx))
  files = files[fileIdx]

# process each file in turn
benchmarks = data.frame()
for (file in files) {
  filename = gsub("\\.txt.*$", "", basename(file))
  message("creating coexpression networks for file ", filename, " (", 
          which(files == file), " of ", length(files), ")")
  tryCatch({
    dat = read.delim(file) %>%
      as.matrix()
    
    # sample a random 1,000 genes
    set.seed(which(files == file))
    n_genes = 1000
    if (ncol(dat) < n_genes)
      next
    dat = dat[, sample(ncol(dat), n_genes)]
    
    # benchmark each correlation
    coefs = dismay::metrics()
    benchmark = data.frame()
    for (coef in coefs) {
      message("  processing coefficient ", coef, " ...")
      mb = microbenchmark(suppressMessages(dismay::dismay(dat, metric = coef)),
                          times = 10)
      benchmark %<>% rbind(data.frame(
        file = filename, coef = coef, mean = mean(mb$time), 
        median = median(mb$time)))
    }
    benchmarks %<>% rbind(benchmark)
  }, error = function(e) {
    message("     ", e)
  })
}

# write benchmark
output_dir = "data/benchmark"
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)
output = ifelse(is.null(fileIdx),
                file.path(output_dir, "benchmark.txt"),
                file.path(output_dir, paste0("benchmark-", fileIdx, ".txt")))
write.table(benchmarks, output, quote = F, sep = "\t", row.names = F)
