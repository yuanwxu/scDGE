log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

name_converter <- function(name_vector) {
    # Helper function to clean cell types names for use in R
    name_vector <- gsub("\\+", ".", name_vector) # replace "+" with "."
    name_vector <- gsub(" ", "_", name_vector) # replace " " with "_"
    name_vector <- gsub("/", "_", name_vector) # replace "/" with "_"
    name_vector <- gsub("-", ".", name_vector)  # replace "-" with "."
    return(name_vector)
}

# Load data
counts <- read.csv(snakemake@input[["counts"]], row.names = 1, check.names = FALSE)
coldata <- read.csv(snakemake@input[["coldata"]], row.names = 1, check.names = FALSE)

# Specify output directory
out_dir <- snakemake@output[["out_dir"]]
dir.create(out_dir, recursive = TRUE)

# Ensure column order matches
stopifnot(all(colnames(counts) == rownames(coldata)))

# Clean up cell type names to avoid potential R naming issues
coldata$cell_type <- name_converter(coldata$cell_type)

# Restrict to a subset of cell types if specified
celltype_pat <- snakemake@params[["celltype_pat"]]

if (!is.null(celltype_pat)) {
    patl <- grepl(celltype_pat, colnames(counts))
    counts <- counts[, patl, drop = FALSE]
    coldata <- coldata[patl, ]
}

# Get unique cell types
cell_types <- unique(coldata$cell_type)

for (ct_idx in seq_along(cell_types)) {
    ct <- cell_types[ct_idx]
    cat("Processing cell type:", ct, "\n")

    # Subset counts and coldata for this cell type
    ct_samples <- rownames(coldata[coldata$cell_type == ct, ])
    counts_ct <- counts[, ct_samples, drop = FALSE]
    coldata_ct <- coldata[ct_samples, ]

    # Check if we have both conditions and enough replicates
    conditions_present <- table(coldata_ct$condition)

    if (length(conditions_present) < 2) {
        cat("Skipping - only one condition present for ", ct, "\n")
        next
    }

    if (any(conditions_present < snakemake@params[["min_reps_per_condition"]])) {
        cat("Skipping - insufficient replicates per condition for ", ct, "\n")
        next
    }

    # Write 
    write.csv(counts_ct, file.path(out_dir, paste0("counts_", ct_idx, ".csv")))
    write.csv(coldata_ct, file.path(out_dir, paste0("coldata_", ct_idx, ".csv")))
}



    