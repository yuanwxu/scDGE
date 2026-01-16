# Capture output and error messages to log file
log <- file(snakemake@log[[1]], open="wt")
sink(log, type="output")
sink(log, type="message")

library(tidyverse)
library(ggrepel)
library(pheatmap)
library(EnhancedVolcano)
library(RColorBrewer)

source("scripts/deseq2_helpers.r")

res_combined <- read_csv(snakemake@input[["de_res"]])
counts_normalized <- read.csv(snakemake@input[["counts_norm"]], row.names = 1, check.names = FALSE)
coldata <- read.csv(snakemake@input[["coldata"]], row.names = 1, check.names = FALSE)

# Ensure column order matches
stopifnot(all(colnames(counts_normalized) == rownames(coldata)))

# Clean up cell type names to avoid potential R naming issues
coldata$cell_type <- name_converter(coldata$cell_type)

# Specify output directory
out_dir <- snakemake@output[["plot_dir"]]
dir.create(out_dir, recursive = TRUE)

# ========== Count DE genes per cell type ==================================
de_counts <- res_combined %>%
  filter(sig != "non-sig") %>%
  group_by(cell_type, sig) %>%
  summarise(count = n(), .groups = "drop")

print("DE gene counts per cell type: ")
print(de_counts %>% tidyr::pivot_wider(names_from = sig, values_from = count, values_fill = 0))


# ========== Bar plot: Number of DE genes per cell type ====================

p1 <- ggplot(de_counts, aes(x = reorder(cell_type, count), y = count, fill = sig)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_fill_manual(values = c("up" = "#E64B35", "down" = "#4DBBD5")) +
  labs(title = "Number of DE genes per cell type",
       x = "Cell type", y = "Number of DE genes", fill = "Direction") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))

ggsave(file.path(out_dir, "bar_degs_by_celltype.pdf"), p1, width = 8, height = 6)


# ========== Volcano plots per cell type =========================

de_counts_total <- de_counts %>%
    group_by(cell_type) %>%
    summarise(total = sum(count))

# Select cell types with at least certain number of DE genes
min_degs <- snakemake@params[["min_degs"]]
top_cell_types <- de_counts_total %>%
    filter(total >= min_degs) %>%
    arrange(desc(total)) %>%
    pull(cell_type)

# Get significance thresholds
padj_level <- snakemake@params[["padj_level"]]
abs_lfc_min <- snakemake@params[["abs_lfc_min"]]

for (ct in top_cell_types) {
  
  res_ct <- res_combined %>% filter(cell_type == ct)
  
  p_volcano <- EnhancedVolcano(
    res_ct,
    lab = res_ct$gene,
    x = 'log2FoldChange',
    y = 'padj',
    title = paste0('Volcano plot: ', ct),
    pCutoff = padj_level,
    FCcutoff = abs_lfc_min,
    pointSize = 2.0,
    labSize = 3.0,
    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
    colAlpha = 0.6,
    legendPosition = 'right',
    legendLabSize = 10,
    legendIconSize = 3.0,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    max.overlaps = 20
  )
  
  ggsave(file.path(out_dir, paste0("volcano_", gsub("[^[:alnum:]]", "_", ct), ".pdf")), 
         p_volcano, width = 10, height = 8)
}


# ======== Heatmaps of top DE genes across cell types =================================

# Get top genes by adjustd p-value
n_genes_max <- snakemake@params[["top_n_genes"]]
top_genes_overall <- res_combined %>%
  filter(sig != "non-sig") %>%
  arrange(padj) %>%
  head(n_genes_max) %>%
  pull(gene) %>%
  unique()

# All cell types with DE genes
cell_types_with_de <- de_counts_total$cell_type

# Prepare pheatmap input
pheatmap_input_list <- prepare_pheatmap(counts_normalized, coldata, top_genes_overall,
                                      which.samples = coldata$cell_type %in% cell_types_with_de)

# Find gap positions (between conditions)
gap_positions <- which(diff(as.numeric(pheatmap_input_list$annotation$condition)) != 0)

# Set reference condition for coloring
ref_level <- snakemake@params[["ref_level"]]
pheatmap_input_list$annotation$condition <- relevel(pheatmap_input_list$annotation$condition, ref = ref_level)
cond_values <- levels(pheatmap_input_list$annotation$condition)
n_col_needed <- length(cond_values)

# Create a palette
col_palette <- brewer.pal(max(n_col_needed, 3), "Pastel2")

# Create heatmap
pheatmap(pheatmap_input_list$counts_subset,
         annotation_col = pheatmap_input_list$annotation,
         annotation_colors = list(condition = setNames(col_palette[1:n_col_needed], cond_values)),
         show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         gaps_col = gap_positions,  # Add visual gap between conditions
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
         fontsize_row = 10,
         main = paste0("Top ", length(top_genes_overall), " DE genes across cell types"),
         filename = file.path(out_dir, "heatmap_top_genes_across_celltypes.pdf"),
         width = 12,
         height = 10)

# ========= Heatmaps of top DE genes per cell type ============================

for (ct in cell_types_with_de) {
  
  # Get significant DE genes for this cell type, sorted by log2FC
  sig_genes_df_ct <- res_combined %>%
    filter(cell_type == ct & sig != "non-sig") %>%
    arrange(log2FoldChange)  # Sort by LFC for better visualization
  
  sig_genes_ct <- sig_genes_df_ct$gene
  
  if (length(sig_genes_ct) == 0) {
    cat("  No significant genes. Skipping.\n")
    next
  }
  
  # Limit to top genes
  if (length(sig_genes_ct) > n_genes_max) {
    sig_genes_ct <- sig_genes_ct[1:n_genes_max]
  }

  # Count genes up vs down
  n_up <- sum(sig_genes_df_ct$sig == "up")
  n_down <- sum(sig_genes_df_ct$sig == "down")

  # Prepare pheatmap input
  pheatmap_input_list <- prepare_pheatmap(counts_normalized, coldata, sig_genes_ct,
                                        which.samples = coldata$cell_type == ct, 
                                        annotate_celltype=FALSE)

  # Find gap positions (between conditions)
  gap_positions <- which(diff(as.numeric(pheatmap_input_list$annotation$condition)) != 0)

  # Create heatmap
  ct_safe <- gsub("[^[:alnum:]]", "_", ct)
  pheatmap(pheatmap_input_list$counts_subset,
          annotation_col = pheatmap_input_list$annotation,
          annotation_colors = list(condition = setNames(col_palette[1:n_col_needed], cond_values)),
          show_colnames = FALSE,
          show_rownames = TRUE,
          cluster_rows = FALSE,  # Don't cluster rows (keep sorted by LFC)
          cluster_cols = FALSE,  # Don't cluster columns (keep sorted by condition)
          gaps_col = gap_positions,  # Add visual gap between conditions
          color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
          fontsize_row = 10,
          main = paste0(ct, ": ", n_down, " down, ", n_up, " up (n=", length(sig_genes_ct), " genes)"),
          filename = file.path(out_dir, paste0("heatmap_top_genes_", ct_safe, ".pdf")),
          width = 10,
          height = max(6, length(sig_genes_ct) / 4)
  )
  cat("Heatmap for cell type:", ct, "\n")
  cat("  Up-regulated:", n_up, ", Down-regulated:", n_down, "\n")
}


# ====== Dotplot of top DE genes across cell types ============================

# Select top x DE genes per cell type
max_degs_per_ct <- snakemake@params[["max_degs"]]
top_genes_per_ct <- res_combined %>%
  filter(sig != "non-sig") %>%
  group_by(cell_type) %>%
  arrange(padj) %>%
  slice_head(n = max_degs_per_ct) %>%
  ungroup()

p_dot <- ggplot(top_genes_per_ct, 
                aes(x = cell_type, y = gene, size = -log10(padj), color = log2FoldChange)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size_continuous(range = c(2, 8)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold")) +
  labs(title = "Top DE genes across cell types",
       x = "Cell type", y = "Gene",
       size = "-log10(padj)", color = "log2FC")

ggsave(file.path(out_dir, "dotplot_top_genes_across_celltypes.pdf"), p_dot, width = 12, height = 10)


# ====== MA plots  =========================================================

for (ct in cell_types_with_de) {
  
  res_ct <- res_combined %>% filter(cell_type == ct)
  
  p_ma <- ggplot(res_ct, aes(x = baseMean, y = log2FoldChange, color = sig)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_x_log10() +
    scale_color_manual(values = c("up" = "red", "down" = "blue", "non-sig" = "grey")) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black") +
    labs(title = paste0("MA plot: ", ct),
         x = "Mean expression (log10)", y = "log2 Fold Change") +
    theme_minimal()

  ggsave(file.path(out_dir, paste0("ma_plot_", gsub("[^[:alnum:]]", "_", ct), ".pdf")), 
         p_ma, width = 8, height = 6)
}