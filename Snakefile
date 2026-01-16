import os

configfile: "config.yaml"

rule all:
    input:
        "output/deseq2/plots/"


rule psuedobulk:
    input:
        config['adata']
    output:
        "output/deseq2/counts.csv",
        "output/deseq2/coldata.csv",
        "output/" + os.path.splitext(os.path.basename(config['adata']))[0] + "_pb.h5ad",
        "output/pseudobulk_eda.pdf" 
    params:
        donor_key=config['keys']['donor'],
        condition_key=config['keys']['condition'],
        celltype_key=config['keys']['celltype'],
        covars=config['additional_covariate'],
        n_cells_min=30,
    conda:
        "envs/scanpy.yaml"
    log:
        "logs/pseudobulk/pseudobulk.log"
    script:
        "scripts/create_pseudobulk.py"


rule normalize_counts:
    input:
        counts="output/deseq2/counts.csv",
        coldata="output/deseq2/coldata.csv"
    output:
        "output/deseq2/counts_normalized.csv"
    params:
        min_samples=4, # keep genes expressed in at least `min_samples` samples with >= `min_counts` counts
        min_counts=10,
    conda:
        "envs/r_deseq2.yaml"
    log:
        "logs/deseq2/normalize_counts.log"
    script:
        "scripts/normalize_counts.R"


checkpoint split_by_celltype:
    input:
        counts="output/deseq2/counts.csv",
        coldata="output/deseq2/coldata.csv"
    output:
        out_dir=temp(directory("output/deseq2/by_celltype/")) # delete at the end
    params:
        celltype_pat=config['deseq2']['subset_celltype_pattern'], # pattern to subset cell types
        min_reps_per_condition=config['deseq2']['min_reps_per_condition'],
    conda:
        "envs/r_deseq2.yaml"
    log:
        "logs/deseq2/split_by_celltype.log"
    script:
        "scripts/split_by_celltype.R"


rule run_deseq2:
    input: 
        counts="output/deseq2/by_celltype/counts_{ct}.csv",
        coldata="output/deseq2/by_celltype/coldata_{ct}.csv"
    output:
        "output/deseq2/by_celltype_results/res_{ct}.csv"
    params:
        ref_level=config['deseq2']['ref_condition'],
        large_n=20,
        lfc_shrink="apeglm",
        padj_level=config['deseq2']['padj_level'],
        abs_lfc_min=config['deseq2']['abs_lfc_min'],
    conda:
        "envs/r_deseq2.yaml"
    log:
        "logs/deseq2/deseq2_{ct}.log"
    threads: 8
    script:
        "scripts/run_deseq2.R"


def get_deseq2_outputs(wildcards):
    # This triggers after the checkpoint
    checkpoint_output = checkpoints.split_by_celltype.get(**wildcards).output[0]
    # Find what files R actually created
    import glob
    cts = [f.split("_")[-1].replace(".csv", "") 
           for f in glob.glob(f"{checkpoint_output}/counts_*.csv")]
    return expand("output/deseq2/by_celltype_results/res_{ct}.csv", ct=cts)

rule gather_deseq2:
    input:
        get_deseq2_outputs
    output:
        protected("output/deseq2/deseq2_results.csv")
    conda:
        "envs/r_deseq2.yaml"
    log:
        "logs/deseq2/gather_deseq2.log"
    script:
        "scripts/gather_deseq2.R"


rule plot_deseq2:
    input:
        de_res="output/deseq2/deseq2_results.csv",
        counts_norm="output/deseq2/counts_normalized.csv",
        coldata="output/deseq2/coldata.csv"
    output:
        plot_dir = directory("output/deseq2/plots/")
    params:
        padj_level=config['deseq2']['padj_level'],
        abs_lfc_min=config['deseq2']['abs_lfc_min'],
        ref_level=config['deseq2']['ref_condition'],
        min_degs=config['deseq2']['plot']['min_degs_per_celltype'],
        max_degs=config['deseq2']['plot']['max_degs_per_celltype'],
        top_n_genes=config['deseq2']['plot']['heatmap_top_n_genes'],
    conda:
        "envs/r_deseq2.yaml"
    log:
        "logs/deseq2/plot_deseq2.log"
    script:
        "scripts/plot_deseq2.R"
