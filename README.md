## Evaluating measures of association for single-cell transcriptomics

This directory contains R code and data required to reproduce and extend the analyses presented in the paper, "Evaluating measures of association for single-cell transcriptomics."

### Data preprocessing

The analysis includes 211 single-cell RNA-seq (scRNA-seq) datasets, of which 162 were obtained from the Gene Expression Omnibus (GEO), 10 were obtained from the 10X Genomics [website](https://support.10xgenomics.com/single-cell-gene-expression/datasets), and 39 were obtained from [mousebrain.org](http://mousebrain.org). 

The raw files obtained from the GEO are located in `data/geo/raw`. The raw files obtained from 10xgenomics.com and mousebrain.org are not included in this repository, due to their size. However, processed and filtered files from the 10X Genomics website are provided in `data/10xgenomics/com/processed` and `data/10xgenomics/com/filtered`, respectively. Only filtered files from mousebrain.org are included, in `data/loom/filtered`. 

GEO files were preprocessed into a common format using the scripts in the `R/geo/preprocessing` folder. These were subsequently filtered to exclude genes that were not detectably expressed in 95% of cells or more, as well as to exclude non-protein-coding genes, using the `R/geo/filter-geo.R` script. The R scripts used to process and filter the other two types of datasets are available in `R/10xgenomics.com` and `R/mousebrain.org`, respectively. 

### Coexpression network generation

Coexpression networks were generated from filtered gene expression files with seventeen different measures of association as implemented in the `dismay` R package. The script `R/coexpr/write-matrices.R` writes coexpression matrices. These were subsequently filtered to exclude genes absent from 80% of cells or more for the main analyses, as well as at different levels of filtering for two supplementary analyses, using the `R/coexpr/filter-matrices.R` script. The `dismay` R package is available from [GitHub](https://github.com/skinnider/dismay).

### Functional coherence analysis

The functional coherence of each network was calculated for each GO term using the `R/function/calculate-auroc.R` script. The outputs were subsequently consolidated (`consolidate-auroc.R` for the main analysis and `consolidate-auroc-filtered.R` for the supplementary analyses), and figures and statistical tests were performed in the `R/function/plot-auroc*` scripts. The contribution of different experimental and analytical factors to functional coherence was assessed with univariate linear models in the `R/function/analyze-r2.R` script. For the most part, our analysis focused on one randomly sampled dataset from each publication except where noted otherwise in the paper, and this random sample was held consistent for all analyses; it is provided in `data/geo/one-dataset-per-publication.txt`. 

### Network overlap analysis

Four different types of biological networks were analyzed: protein-protein interactions from [HIPPIE](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/information.php), signalling networks from [OmniPath](http://omnipathdb.org/), metabolic pathway co-membership networks from [Reactome](https://reactome.org/), and gene co-occurrence networks derived from text mining from [STRING](https://string-db.org/). All four are provided in their respective directories within `data/networks` with preprocessing code provided in `R/networks/rewire-*`. To provide a null model for the observed overlap, each network was rewired 1,000 times using the `R/networks/rewire-networks.R` script; because of the size of these files, they are not included in the repository. The observed overlap was compared to random expectation using the script `R/networks/calculate-overlap.R`, and visualized using `R/networks/plot-network-overlap.R`. 

### Cell clustering analysis

Cell clustering with each of the 17 measures of association was analyzed for the [Li et al.](https://www.nature.com/articles/ng.3818) dataset, using two different algorithms (hierarchical clustering and Louvain clustering of the shared-nearest-neighbor graph). The cell-cell matrices was generated with `R/clustering/write-cell-matrices.R` and the principal analysis was performed using `R/clustering/analyze-hclust.R` and `R/clustering/analyze-snn.R`, respectively. The dendrograms obtained by hierarchical clustering with each method were also plotted with both the cell line of origin and batch in the script `R/clustering/plot-dendrograms.R`. 

### Reproducibility analysis

To analyze the reproducibility of coexpression networks inferred with each measure of association we made use of five scRNA-seq datasets of human pancreatic alpha, beta, and delta cells, (all of which were obtained from the GEO and are in `data/geo/filtered` directory). We considered all 30 combinations of each cell type across five datasets. For each pair of datasets, we filtered both coexpression matrices to the intersect of the genes present in either network, then calculated the Spearman correlation between the two matrices. To assess statistical significance, we subsequently permuted the matrices 100 times, following the permutation procedure of the Mantel test. We then obtained the z score of the observed Spearman correlation relative to random expectation, and visualized the results in `R/reproducibility/plot-reproducibility.R`. 

### Disease gene analysis

We analyzed disease gene prediction for CNS-related disorders and cell type-specific coexpression of genes associated with cerebrovascular disease; both lists of disease genes were obtained from [Phenopedia](https://phgkb.cdc.gov/PHGKB/startPagePhenoPedia.action), and preprocessed using the scripts `R/disease/preprocess-phenopedia.R` and `R/disease/create-mesh-map.R`. We subsequently used the same framework as in the functional coherence analysis (reimplemented in `R/disease/analyze-auroc-disease.R`) to analyze disease gene coexpression in each case, visualizing the results with the script `R/disease/plot-auroc-disease.R`. 
## SCT-MoA
