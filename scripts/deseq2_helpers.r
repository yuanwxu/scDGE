# Define helper functions for DESeq2 analysis


name_converter <- function(name_vector) {
    # Helper function to clean cell types names for use in R
    name_vector <- gsub("\\+", ".", name_vector) # replace "+" with "."
    name_vector <- gsub(" ", "_", name_vector) # replace " " with "_"
    name_vector <- gsub("/", "_", name_vector) # replace "/" with "_"
    name_vector <- gsub("-", ".", name_vector)  # replace "-" with "."
    return(name_vector)
}


filter_counts <- function(counts, coldata, min_samples=4, min_counts=10, min_total_counts=NULL) {
    # Filter raw count matrix

    # Filter genes: keep genes expressed in at least `min_samples` samples 
    # with >= `min_counts` counts
    genes_to_keep <- rowSums(counts >= min_counts) >= min_samples
    counts_filtered <- counts[genes_to_keep, ]

    # Optionally filter samples with low total counts
    if (!is.null(min_total_counts)) {
        samples_to_keep <- colSums(counts_filtered) >= min_total_counts
        counts_filtered <- counts_filtered[, samples_to_keep]
        coldata_filtered <- coldata[samples_to_keep, ]

        return(list(counts=counts_filtered, coldata=coldata_filtered))
    }
    
    return(list(counts=counts_filtered, coldata=coldata))
}


normalize_counts <- function(counts, coldata, normalization="vst", design = ~ condition + cell_type) {
  # Normalize raw count data
  # 
  # design: the design formula to construct a DESeqDataSet object and seen by vst
  if (normalization == "vst") { 
    # Use 'VST' transformation in DESeq2, 'rlog' is too slow
    all_dds <- DESeqDataSetFromMatrix(
      countData = counts,
      colData = coldata,
      design = design
      )
    all_dds <- estimateSizeFactors(all_dds)
    vst_counts <- vst(all_dds, blind = FALSE)
    counts_norm <- assay(vst_counts)
  } 
  else {
    stop("Not implemented, only 'vst' normalization is supported.")
  }

  return (counts_norm)
} 


prepare_pheatmap <- function(counts_norm, coldata, de_genes, which.samples, sort_by_condition=TRUE,
                            annotate_celltype=TRUE) {
  # Get input to pheatmap for plotting top DE genes, expects normalized counts
  # which.samples: a logical vector of which samples to show in the heatmap

  if (!any(which.samples)) {
    stop("No samples selected for heatmap.")
  }

  # Extract normalized counts for the top genes
  if (sort_by_condition) {
    # Sort columns by condition for visual separation
    condition_vec <- coldata$condition[which.samples]
    sort_idx <- order(condition_vec)

    counts_subset <- counts_norm[de_genes, which.samples, drop=FALSE][, sort_idx, drop=FALSE]

    annotation_col <- data.frame(
      condition = factor(condition_vec[sort_idx]),
      row.names = colnames(counts_subset)
    )

    if (annotate_celltype) {
        annotation_col["cell_type"] <- coldata$cell_type[which.samples][sort_idx]
    }

  } else {
    counts_subset <- counts_norm[de_genes, which.samples, drop=FALSE]

    annotation_col <- data.frame(
      condition = factor(coldata$condition[which.samples]),
      row.names = colnames(counts_subset)
    )

    if (annotate_celltype) {
        annotation_col["cell_type"] <- coldata$cell_type[which.samples]
    }
  }

  return(list(counts_subset=counts_subset, annotation=annotation_col))
}


prepare_deseq <- function(counts, coldata, ref_level, min_reps_per_condition=2,
                        include_design_donor_id=FALSE, ...){
  # Basic QC and filtering for DESeq2, for a single cell type
  #
  # Parameters:
  # counts: count matrix (genes x samples)
  # coldata: data frame with sample metadata
  # min_reps_per_condition: minimum replicates per condition 
  # include_design_donor_id: whether to include donor/sample ID in design matrix,
  #                          works only if there are multiple biological replicates
  #                          per donor, or same donor with different conditions 
  #                          (e.g., paired/before-after samples). If donor and
  #                          condition are perfectly confounded, this will cause
  #                          the model matrix to be singular.
  # ... additional arguments to edgeR::filterByExpr

  # Ensure column order matches
  stopifnot(all(colnames(counts) == rownames(coldata)))

  # Check that only one cell type is present
  if (length(unique(coldata$cell_type)) != 1) {
    stop("For fitting DESeq model separately for each cell type, expected only 
      one cell type in coldata. Found: ", unique(coldata$cell_type))
  }

  # Check if we have both conditions and enough replicates
  conditions_present <- table(coldata$condition)

  if (any(conditions_present < min_reps_per_condition)) {
    stop("Insufficient replicates per condition for cell type ",
         coldata$cell_type[1], ": ", 
         paste(names(conditions_present), conditions_present, sep="=", collapse=", "))
  }

  # Set reference level
  coldata$condition <- relevel(factor(coldata$condition), ref = ref_level)

  # Filter low-count genes
  if (include_design_donor_id) {
    des_mat <- model.matrix(~ condition + sample, data = coldata)
  } else {
    des_mat <- model.matrix(~ condition, data = coldata)
  }
  
  # Use edgeR::filterByExpr logic for filtering
  keepl <- edgeR::filterByExpr(counts, design = des_mat, ...)
  counts_filtered <- counts[keepl, ]
  cat("Filtered from", nrow(counts), "to", nrow(counts_filtered), "genes after low-count filtering\n")

  # Create DESeqDataSet
  if (include_design_donor_id) {
    dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                                  colData = coldata,
                                  design = ~ sample + condition)
  } else {
    dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                                  colData = coldata,
                                  design = ~ condition)
  }
  return(dds) 
}


prepare_deseq_interaction <- function(counts_pb, coldata, min_reps_per_condition=2, ...){
    # Basic QC and filtering for DESeq2, using an interaction model (cell_type:condition)
    #
    # Parameters:
    # counts_pb: pseudobulk count matrix (genes x samples)
    # coldata: data frame with sample metadata
    # min_reps_per_condition: minimum replicates per condition
    # ... additional arguments to edgeR::filterByExpr

    # Ensure column order matches
    stopifnot(all(colnames(counts) == rownames(coldata)))

    # Check number of samples within each cell_type × condition combinations 
    design_table <- table(coldata$cell_type, coldata$condition)
    # print(design_table)
    
    # Filter out cell types with insufficient samples (R or NR)
    keepl_ct <- apply(design_table, 1, function(x) all(x >= min_reps_per_condition))
    coldata_filtered <- coldata %>% 
        filter(cell_type %in% names(keepl_ct)[keepl_ct])
    cat("Filtered from", nrow(coldata), "to", nrow(coldata_filtered), 
        "samples/pseudobulk groups after removing cell types with insufficient replicates\n")
    
    counts_filtered <- counts_pb[, rownames(coldata_filtered)]

    # Filter low-count genes
    des_mat <- model.matrix(~ cell_type + condition + cell_type:condition, data=coldata_filtered)

    # Use edgeR::filterByExpr logic for filtering
    keepl <- edgeR::filterByExpr(counts_filtered, design = des_mat, ...)
    counts_filtered <- counts_filtered[keepl, ]
    cat("Filtered from", nrow(counts_pb), "to", nrow(counts_filtered), "genes after low-count filtering\n")

    # Create DESeqDataset with interaction design
    dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                                    colData = coldata_filtered,
                                    design = ~ cell_type + condition + cell_type:condition)
    
    return(dds)
}


build_contrast_specs <- function(cell_types, ref_cell_type, result_names) {
    specs <- list()
    # Function to build contrast specs programmatically in an interaction model
    # Return a list of contrast specifications for each cell type (main condition effect + interaction)

    # For reference cell type, look for the main condition effect
    # Pattern: "condition_NR_vs_R" or similar
    cond_names <- grep("^condition_", result_names, value = TRUE)
    if (length(cond_names) != 1)
        stop("Expected exactly one main condition effect, found:", length(cond_names))
  
    for (ct in cell_types) {
        if (ct == ref_cell_type) {
            specs[[ct]] <- list(type = "name", value = cond_names)
            cat("Cell type:", ct, "(reference) -> using:", cond_names, "\n")

        } else {
        # For non-reference cell types, look for interaction terms
        # Pattern: "cell_type*interaction_term" or similar
        
        # Try to find a coefficient that matches this cell type and interaction
        interaction_pattern <- paste0(".*", ct, ".condition.*")
        interaction_names <- grep(interaction_pattern, result_names, value = TRUE)

        if (length(interaction_names) != 1) {
            stop("Expected exactly one interaction term for cell type ", ct, 
                 ", found: ", length(interaction_names), "(", paste(interaction_names, collapse = ", "), ")")
        }
        
        # Use numeric contrast by finding indices
        # interaction_idx <- match(interaction_names, result_names)
        
        # Build a contrast vector with 1 for this interaction
        # contrast_vec <- rep(0, length(result_names))
        # contrast_vec[interaction_idx] <- 1
        
        # Also add the main condition effect
        # contrast_vec[match(cond_names, result_names)] <- 1
        
        # specs[[ct]] <- list(type = "contrast", value = contrast_vec)

        # Build contrast list: main condition effect + interaction effect
        contrast_list <- list(c(interaction_names, cond_names))
        specs[[ct]] <- list(type = "contrast", value = contrast_list)
        }
    }
    return(specs)
}


extract_deseq_results_interaction <- function(dds, cell_types, contrast_specs, lfc_shrink=TRUE) {
    # Function to extract DESeq2 results per cell type in an interaction model
    # Return results of all cell types in a data frame
    res_list <- list()

    for (ct in cell_types) {
        spec <- contrast_specs[[ct]]
    
        # Get results based on type (either by name or by contrast vector)
        if (spec$type == "name") { # ref cell type
            res <- results(dds, name = spec$value)
            if (lfc_shrink) {
                # Apply shrinkage 
                res <- lfcShrink(dds, coef = spec$value, res = res, type = "ashr")
            }
        } else if (spec$type == "contrast") {
            res <- results(dds, contrast = spec$value)
            if (lfc_shrink) {
                # Apply shrinkage (type='apeglm' does not support "contrast")
                res <- lfcShrink(dds, contrast = spec$value, res = res, type = "ashr")
            }
        }
        
        # Convert to data frame and add cell type info
        res_df <- as.data.frame(res) %>%
            rownames_to_column("gene") %>%
            mutate(
                cell_type = ct,
                sig = case_when(
                !is.na(padj) & padj < padj_level & log2FoldChange > abs_lfc_min ~ "up",
                !is.na(padj) & padj < padj_level & log2FoldChange < -abs_lfc_min ~ "down",
                TRUE ~ "non-sig"
                )
            )
        
        res_list[[ct]] <- res_df
        
        # Print summary
        cat("DE summary for", ct, ":\n")
        cat("Up-regulated (padj < ", padj_level, ", log2FC >", abs_lfc_min, "):", sum(res_df$sig == "up"), "\n")
        cat("Down-regulated (padj < ", padj_level, ", log2FC < ", -abs_lfc_min, "):", sum(res_df$sig == "down"), "\n")
        cat("Non-significant:", sum(res_df$sig == "non-sig"), "\n\n")
    }

    # Combine all results into a single data frame
    res_combined <- bind_rows(res_list)

    return(res_combined)
}