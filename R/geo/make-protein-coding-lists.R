# Create lists of protein-coding genes in human and mouse. 
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")
library(tidyverse)
library(magrittr)
library(rtracklayer)

# read human and mouse GTF files 
hsa = as.data.frame(import("data/ensembl/Homo_sapiens.GRCh38.91.gtf.gz"))
mmu = as.data.frame(import("data/ensembl/Mus_musculus.GRCm38.91.gtf.gz"))

# filter to protein-coding genes from both
hsa_genes = hsa %>%
  dplyr::filter(gene_biotype == 'protein_coding' & type == 'gene') %>%
  pull(gene_id) %>%
  unique()
mmu_genes = mmu %>%
  dplyr::filter(gene_biotype == 'protein_coding' & type == 'gene') %>%
  pull(gene_id) %>%
  unique()

# combine
dat = rbind(data.frame(source = 'human', gene = hsa_genes),
            data.frame(source = 'mouse', gene = mmu_genes))

# write
write_and_gzip(dat, "data/ensembl/protein_coding_genes.txt")
