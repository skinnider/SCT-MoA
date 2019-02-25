# Dataset: GSE67835
# Cells: single cells from the human CNS
# Values: counts
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")
library(GEOquery)

# define accession
accession = "GSE67835"

# read SOFT file
soft = getGEO(filename = "data/geo/raw/GSE67835/GSE67835_family.soft")
gsms = GSMList(soft)
header = gsms[[1]]@header

# format header as data frame
characteristics = header$characteristics_ch1
header = data.frame(header[!names(header) %in% c(
  "characteristics_ch1", "contact_country", "data_row_count", 
  "data_processing", "extract_protocol_ch1", "relation")])
cellTypes = characteristics[grepl("cell type", characteristics)]
cellTypes = gsub("cell type: ", "", cellTypes)
header$cellType = cellTypes

# ignore fetal cells and types with < 30 cells
cellTypeCounts = table(cellTypes)
validCellTypes = names(cellTypeCounts)[cellTypeCounts >= 30]
validCellTypes = validCellTypes[!grepl("fetal", validCellTypes)]
validCellTypes = validCellTypes[validCellTypes != "hybrid"] ## ignore hybrid
header = header[header$cellType %in% validCellTypes,]

# read files from each cell type in turn
header$file = basename(header$supplementary_file_1)
for (cellType in validCellTypes) {
  message(".. preprocessing ", cellType, " ...")
  # subset the header
  subset = header[header$cellType == cellType,]
  # read data
  files = file.path("data/geo/raw", accession, subset$file)
  dat = lapply(files, function(file) {
    # read file
    tmp = read.delim(file, header = F)
    # trim whitespace
    tmp[, 1] = trimws(tmp[, 1])
    tmp
  })
  dat = do.call(data.frame, dat)
  genes = dat[, 1]
  # remove unneeded columns
  dat = dat[, grepl("V2", colnames(dat))]
  
  # map to ensembl
  dat0 = map_genes(dat, genes, from = "ALIAS", to = "ENSEMBL",
                   db = org.Hs.eg.db)
  
  # convert to matrix and transpose
  dat1 = t(dat0[, -1])
  colnames(dat1) = dat0[,1]
  rownames(dat1) = gsub("_.*$", "", basename(files))
  
  # write
  message("... writing ", cellType, " with ", nrow(dat1), " cells ...")
  output = paste0("data/geo/processed/", accession, "_", cellType, ".txt")
  write_and_gzip(dat1, output)
}
