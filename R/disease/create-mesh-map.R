# Create map of CUI to MESH terms from the Human Disease Ontology
# DOID.csv.gz download at: https://bioportal.bioontology.org/ontologies/DOID/
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
source("R/functions.R")
library(tidyverse)

# Get MeSH term - CUI term map
map = read.csv("data/disease/phenopedia/DOID.csv.gz") %>%
  mutate(CUI = gsub(".*\\|UMLS_CUI:", "", database_cross_reference)) %>%
  mutate(CUI = gsub("\\|.*", "", CUI)) %>%
  mutate(CUI = ifelse(grepl(":", CUI), NA, CUI)) %>%
  mutate(mesh = gsub(".*\\|MESH:", "", database_cross_reference)) %>%
  mutate(mesh = gsub("\\|.*", "", mesh)) %>%
  mutate(mesh = ifelse(grepl(":", mesh), NA, mesh)) %>%
  dplyr::select(CUI, mesh) %>%
  filter(!is.na(CUI)) %>%
  filter(!is.na(mesh)) %>%
  filter(CUI != "") %>%
  filter(mesh != "")
write_and_gzip(map, "data/disease/phenopedia/CUI-MeSH-map.txt")
