# General-purpose functions

map_genes = function(dat, keys, from = "SYMBOL", to = "ENSEMBL", db = NULL) {
  if (is.null(db) | !"OrgDb" %in% class(db))
    stop("argument `db` must be an object of class OrgDb")
  if (is.null(from) | !from %in% keytypes(db))
    stop("invalid `from` argument: ", from)
  if (is.null(to) | !to %in% keytypes(db))
    stop("invalid `to` argument: ", to)
  map = suppressMessages(AnnotationDbi::select(db, keys = keys, columns = to, 
                                               keytype = from))
  values = map[[to]][match(keys, map[[from]])]
  aggregate(dat, by = list(values), FUN = sum)
}

write_and_gzip = function(data, filepath) {
  write.table(data, filepath, sep = "\t", quote = F, row.names = F)
  system(paste("gzip --force", filepath))
}

clean_metric = function(vec) {
  as.character(fct_recode(vec,
                          "Manhattan distance" = "manhattan",
                          "Euclidean distance" = "euclidean",
                          "Canberra distance" = "canberra",
                          "Zero-inflated Kendall correlation" = "zi_kendall",
                          "Kendall correlation" = "kendall",
                          "Spearman correlation" = "spearman",
                          "Pearson correlation" = "pearson",
                          "Jaccard index" = "jaccard",
                          "Dice coefficient" = "dice",
                          "Hamming distance" = "hamming",
                          "Co-dependency index" = "binomial",
                          "Biweight midcorrelation" = "bicor",
                          "Cosine distance" = "cosine",
                          "Weighted rank correlation" = "weighted_rank",
                          "ϕs" = "phi_s",
                          "ρp" = "rho_p"
  ))
}

get_children = function(term) {
  library(MeSH.PCR.db)
  # get term and children
  terms = select(MeSH.PCR.db, keys = term, columns = c("PARENT", "CHILD"), 
                 keytype = "PARENT")
  # get CUI-MeSH map
  map = read.delim("data/disease/phenopedia/CUI-MeSH-map.txt.gz")
  # subset
  sub = map %>% filter(mesh %in% terms$CHILD)
  # return children
  return(unique(sub$CUI))
}
