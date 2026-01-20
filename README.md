# A Differential Gene Expression (DGE) analysis workflow for single cell data

A streamlined, Snakemake-powered pipeline for Differential Gene Expression (DGE) analysis of single cell data using **DESeq2**. This workflow takes an [anndata](https://anndata.readthedocs.io/en/latest/index.html) object containing the cell-by-gene count matrix, cell type labels, and the phenotypic conditions as minimum. It creates pseudobulks and then uses DESeq2 to run per cell type DGE analysis. It generates the new psudobulk anndata object, diagnostic plots for identifying covariates that should be included in the subsequent DESeq2 design matrix, the DESeq2 results of significant genes across cell types, and various plots for further analysis.  


## Installation

### Clone this repository

```bash  
git clone https://github.com/yuanwxu/scDGE.git
cd scDGE
```

### Create a master environment to run the pipeline

Run the setup script:

```bash
./setup.sh
```

## Usage

### Activate the master environment

```bash
conda activate snakemake
```

### Edit the configuration file to suit your data

Find the `config.yaml` in the folder and specify your data, covariates and parameters for analysis.

### Run the pipeline

To preview what will happen (dry run):
```bash
snakemake -np
```

To run the full analysis:
```bash
snakemake --use-conda --cores 8
```

### In case there is an error

Check the `logs/` folder for detailed error message.


## Results

All resutls are saved in the `output/` folder.

- `<mydata>_pb.h5ad`: the pseudobulk anndata object created from the original single cell anndata object `<mydata>.h5ad`.
- `pseudobulk_eda.pdf`: PCA plots of the pseudolbulk data colored by covariates.
- `deseq2/`: this folder contains the main deseq2 results:
    - `counts.csv`, `coldata.csv`: count and metadata input to DESeq2
    - `counts_normalized.csv`: normalized count data for plotting heatmaps
    - **`deseq2_results.csv`**: the DESeq2 results of significant genes across cell types.
    - **`plots/`**: all plots are here.