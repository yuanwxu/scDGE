# Capture output and error messages to log file
log <- file(snakemake@log[[1]], open="wt")
sink(log, type="output")
sink(log, type="message")

library(dplyr)
library(DESeq2)

source("scripts/deseq2_helpers.r")

# Load data
counts <- read.csv(snakemake@input[["counts"]], row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)

coldata <- read.csv(snakemake@input[["coldata"]], row.names = 1, check.names = FALSE)

ct <- coldata$cell_type[1]

# Prepare DESeq2 dataset
dds <- prepare_deseq(counts, coldata, 
                    ref_level=snakemake@params[["ref_level"]], 
                    large.n=snakemake@params[["large_n"]])

# Run DESeq2
if (!is.null(dds)) {
    dds <- DESeq(dds, parallel = TRUE)

    # Extract results
    cond_levels <- levels(colData(dds)$condition)
    res <- results(dds, contrast = c("condition", cond_levels[2], cond_levels[1]))

    # Apply log fold change shrinkage
    coef_names <- resultsNames(dds)
    res <- lfcShrink(dds, 
                    coef = coef_names[which(coef_names == paste0("condition_", cond_levels[2], "_vs_", cond_levels[1]))],
                    res = res, 
                    type = snakemake@params[["lfc_shrink"]])
    
    # Assign significance based on predefined thresholds
    padj_level <- snakemake@params[["padj_level"]]
    abs_lfc_min <- snakemake@params[["abs_lfc_min"]]

    res_df <- as.data.frame(res) %>%
        tibble::rownames_to_column("gene") %>%
        mutate(
            cell_type = ct,
            sig = case_when(
                !is.na(padj) & padj < padj_level & log2FoldChange > abs_lfc_min ~ "up",
                !is.na(padj) & padj < padj_level & log2FoldChange < -abs_lfc_min ~ "down",
                TRUE ~ "non-sig"
            )
        )

    # Print summary
    cat("DE summary for", ct, ":\n")
    cat("Up-regulated (padj < ", padj_level, ", log2FC >", abs_lfc_min, "):", sum(res_df$sig == "up"), "\n")
    cat("Down-regulated (padj < ", padj_level, ", log2FC < ", -abs_lfc_min, "):", sum(res_df$sig == "down"), "\n")
    cat("Non-significant:", sum(res_df$sig == "non-sig"), "\n\n")

    # Write results
    readr::write_csv(res_df, snakemake@output[[1]])
    
} else { # Skipped due to insufficient data
    stop("Skipping cell type:", ct, "due to insufficient data\n\n")
}




