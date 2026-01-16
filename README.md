# A Differential Gene Expression (DGE) analysis workflow for single cell data

A streamlined, Snakemake-powered pipeline for Differential Gene Expression (DGE) analysis of single cell data using **DESeq2**. This workflow takes an [anndata](https://anndata.readthedocs.io/en/latest/index.html) object containing the cell-by-gene count matrix, cell type labels, and the phenotypic conditions as minimum. It creates pseudobulks and then uses DESeq2 to run per cell type DGE analysis. It generates the new psudobulk anndata object, diagnostic plots for identifying covariates that should be included in the subsequent DESeq2 design matrix, the DESeq2 results of significant genes across cell types, and various plots for further analysis.  


## Installation


## Usage