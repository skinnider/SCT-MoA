# Generate 1,000 randomized versions of human and mouse PPI, signalling,
# metabolic, and co-occurrence networks. 
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(igraph)
library(tidyverse)
library(magrittr)
source("R/functions.R")

# usage: rewire-networks.R <idx> <output_dir>
args = commandArgs(trailingOnly = T)
paste(args, collapse = " ")
if (length(args) < 2)
  stop("must provide index of network to rewire [1-8] and output directory")
idx = as.integer(args[1])
output_dir = args[2]
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)

# define all inputs
files = list.files("data/networks", pattern = "*.gz", recursive = T,
                   full.names = T)
files = files[grepl("human|mouse", files)]
dbs = basename(dirname(files))
species = gsub("\\..*$", "", basename(files))
inputs = data.frame(file = files, network = dbs, species = species)

# analyze the network of interest 
network = inputs$network[idx]
species = inputs$species[idx]
file = inputs$file[idx]
set.seed(idx)

# read network
net = drop_na(read.delim(file))
g = graph_from_data_frame(net, directed = F)
m = length(E(g))

# rewire network
bootstraps = 1e3
for (b in seq_len(bootstraps)) {
  if (b %% 10 == 0) 
    message("bootstrap ", b, " of ", bootstraps, " ...")
  rewired = rewire(g, with = keeping_degseq(
    loops = F, niter = 10 * m)) %>%
    igraph::as_data_frame()
  output_file = paste0(network, "-", species, "-", b, ".txt")
  output_path = file.path(output_dir, output_file)
  write_and_gzip(rewired, output_path)
}
