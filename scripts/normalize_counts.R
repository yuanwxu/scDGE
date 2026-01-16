log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(DESeq2)

source("scripts/deseq2_helpers.r")

name_converter <- function(name_vector) {
    # Helper function to clean cell types names for use in R
    name_vector <- gsub("\\+", ".", name_vector) # replace "+" with "."
    name_vector <- gsub(" ", "_", name_vector) # replace " " with "_"
    name_vector <- gsub("/", "_", name_vector) # replace "/" with "_"
    name_vector <- gsub("-", ".", name_vector)  # replace "-" with "."
    return(name_vector)
}

# Load pseudobulk count data
counts <- read.csv(snakemake@input[["counts"]], row.names = 1, check.names = FALSE)
coldata <- read.csv(snakemake@input[["coldata"]], row.names = 1, check.names = FALSE)

# Ensure column order matches
stopifnot(all(colnames(counts) == rownames(coldata)))

# Clean up cell type names to avoid potential R naming issues
coldata$cell_type <- name_converter(coldata$cell_type) 

# Basic filtering for the raw count data before normalization
filter_res <- filter_counts(counts, coldata, 
                            min_samples=snakemake@params[["min_samples"]], 
                            min_counts=snakemake@params[["min_counts"]])
counts_pbmc_filtered <- filter_res$counts
coldata_pbmc_filtered <- filter_res$coldata

# Get normalized counts using VST
counts_norm <- normalize_counts(filter_res$counts, filter_res$coldata)

# Write output
write.csv(counts_norm, snakemake@output[[1]])