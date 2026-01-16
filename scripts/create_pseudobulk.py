import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import scanpy as sc

# DEG analysis per cell type
def create_pseudobulk(adata, donor_key='sample', condition_key='group', celltype_key='celltype3',
                      n_cells_min=30, obs_cols_to_keep=None,
                      counts_deseq2_file="counts_pseudobulk_DESeq2.csv",
                      coldata_deseq2_file="coldata_pseudobulk_DESeq2.csv"):
    """ Create pseudobulk data for DEG analysis. 
    Parameters
    ----------
    adata: The input AnnData object
    donor_key: The key in adata.obs that indicates donor/sample labels
    condition_key: The key in adata.obs that indicates condition labels (e.g., disease vs. control)
    celltype_key: The key in adata.obs that indicates cell type labels
    n_cells_min: Minimum number of cells required per (celltype, sample) group to include in pseudobulk
    obs_cols_to_keep: List of additional columns in adata.obs to keep in the pseudobulk AnnData object,
            should be sample-level variables because the aggregation is done by taking the first value per
            (celltype, sample) group.
    counts_deseq2_file: Output file name for pseudobulk counts matrix for DESeq2
    coldata_deseq2_file: Output file name for pseudobulk colData for DESeq2
    Returns
    -------
    adata_pb: The pseudobulk AnnData object
    """
    # Choose raw counts source
    if adata.raw is not None and adata.raw.X is not None:
        counts = adata.raw.X  # shape: (n_cells, n_genes)
        gene_names = adata.raw.var_names
    elif adata.layers.get('counts') is not None:
        counts = adata.layers['counts']  # shape: (n_cells, n_genes)
        gene_names = adata.var_names
    else:
        raise ValueError("No raw counts available for pseudobulk (need adata.raw or adata.layers['counts']).")
    
    cell_type = adata.obs[celltype_key].astype(str)
    sample = adata.obs[donor_key].astype(str)
    ct_samp = cell_type + '|' + sample
    
    # Keep only pesudobulk groups with at least n_cells_min cells
    obs_df = adata.obs.copy()
    obs_df['ct_samp'] = ct_samp
    group_counts = obs_df.groupby('ct_samp').size()
    groups_to_keep = group_counts[group_counts >= n_cells_min].index.tolist()
    groups_to_drop = [g for g in group_counts.index if g not in groups_to_keep]
    for g in groups_to_drop:
        print(f"Pseudobulk group {g} has less than {n_cells_min} cells ({group_counts[g]} cells). Dropping group {g}.")
    
    if len(groups_to_keep) == 0:
        raise ValueError("No (celltype, sample) groups meet n_cells_min. Adjust n_cells_min.")
    
    # Build pseudobulk counts (genes x pseudobulk groups) for kept groups only
    pb_cols = groups_to_keep  # e.g., "Celltype X|sample1"
    
    if sp.issparse(counts):
        cols = []
        for g in groups_to_keep:
            idx = (ct_samp.values == g)
            col = counts[idx, :].sum(axis=0)  # 1 x n_genes
            col = sp.csr_matrix(col).T # convert and transpose to (n_genes, 1)
            cols.append(col)
        counts_pb = sp.hstack(cols)  # shape: (n_genes, n_keep)
    else:
        cols = []
        for g in groups_to_keep:
            idx = (ct_samp.values == g)
            col = counts[idx, :].sum(axis=0)  # 1 x n_genes
            col = col.reshape(-1, 1)  # transpose to (n_genes, 1)
            cols.append(col)
        counts_pb = np.column_stack(cols)  # shape: (n_genes, n_keep)
        
    # Write counts for DESeq2
    genes = np.asarray(gene_names)
    counts_pb_df = pd.DataFrame(
        counts_pb.toarray() if sp.issparse(counts_pb) else counts_pb,
        index=genes,
        columns=pb_cols
    )
    counts_pb_df.to_csv(counts_deseq2_file)

    # Build colData for kept groups
    coldata = []
    for g in groups_to_keep:
        ct, samp = g.split('|', 1)
        # Get condition for this (celltype, sample) pair
        cond = adata.obs.loc[(adata.obs[celltype_key] == ct) & (adata.obs[donor_key] == samp), condition_key].iloc[0]
        coldata.append({'cell_type': ct, 'sample': samp, 'condition': cond})
    
    coldata_df = pd.DataFrame(coldata, index=groups_to_keep)
    coldata_df.to_csv(coldata_deseq2_file)
    
    print(f"Wrote counts: {counts_deseq2_file}")
    print(f"Wrote colData: {coldata_deseq2_file}")
    
    # Return pseudobulk anndata object
    coldata_df.rename(columns={'cell_type': celltype_key, 
                               'sample': donor_key, 
                               'condition': condition_key}, 
                    inplace=True) # keep the original column names

    if obs_cols_to_keep is not None:
        assert isinstance(obs_cols_to_keep, list), "obs_cols_to_keep must be a list of column names in adata.obs"
        agg_dict = {col: 'first' for col in obs_cols_to_keep if col in adata.obs.columns}
        agg_df = obs_df.groupby('ct_samp').agg(agg_dict)
        agg_df = agg_df.loc[groups_to_keep]
        coldata_df = coldata_df.join(agg_df)
        
    adata_pb = sc.AnnData(counts_pb.T, obs=coldata_df)
    return adata_pb


def pseudobulk_EDA(adata_pb, condition_key, covars, fig_file='pseudobulk_eda.pdf'):
    """ Simple EDA on pseudobulk data to check for covariates possibly confounding
        GE and so should be included in the design matrix in addition to the condition
        of interst for DEG analysis.
        Follows https://www.sc-best-practices.org/conditions/differential_gene_expression.html#pseudobulk
    """
    # Save the raw counts in layer
    adata_pb.layers['counts'] = adata_pb.X.copy()
    
    # Normalize counts and perform PCA
    sc.pp.normalize_total(adata_pb, target_sum=1e6)
    sc.pp.log1p(adata_pb)
    sc.pp.pca(adata_pb)
    
    # Add library size to obs
    adata_pb.obs["lib_size"] = np.sum(adata_pb.layers["counts"], axis=1)
    adata_pb.obs["log_lib_size"] = np.log(adata_pb.obs["lib_size"])
    
    # Visualize PCA colored by different covariates
    sc.pl.pca(adata_pb, color=[condition_key] + covars + ["lib_size", "log_lib_size"], 
              ncols=1, show=False)
    plt.savefig(fig_file, bbox_inches='tight', dpi=300)
    
    # Restore original counts
    adata_pb.X = adata_pb.layers['counts'].copy()


# Create pseudobulk from scRNA-seq for DGE analysis with DESeq2
adata = sc.read_h5ad(snakemake.input[0])
adata_pb = create_pseudobulk(adata, 
                             donor_key=snakemake.params.donor_key,
                             condition_key=snakemake.params.condition_key,
                             celltype_key=snakemake.params.celltype_key,
                             n_cells_min=snakemake.params.n_cells_min,
                             obs_cols_to_keep=snakemake.params.covars,
                             counts_deseq2_file=snakemake.output[0],
                             coldata_deseq2_file=snakemake.output[1])
                            
pseudobulk_EDA(adata_pb, 
               condition_key=snakemake.params.condition_key, 
               covars=snakemake.params.covars, 
               fig_file=snakemake.output[3])

# Save the pseudobulk object
adata_pb.write_h5ad(snakemake.output[2])

