# Load data

dir <- "output/deseq2"

counts <- read.csv(file.path(dir, "counts.csv"), row.names = 1, check.names = FALSE)
coldata <- read.csv(file.path(dir, "coldata.csv"), row.names = 1, check.names = FALSE)

# Ensure column order matches
stopifnot(all(colnames(counts) == rownames(coldata)))

# Split into Tissu and PBMC cell types to run DESeq2 separately
pbmc_pat <- "^PBMC "
pbmcl <- grepl(pbmc_pat, colnames(counts))
counts_pbmc <- counts[, pbmcl]
coldata_pbmc <- coldata[pbmcl, ]

counts_tissue <- counts[, !pbmcl]
coldata_tissue <- coldata[!pbmcl, ]

# Write outputs
dir.create(file.path(dir, "pbmc"), recursive=TRUE)
dir.create(file.path(dir, "tissue"), recursive=TRUE)

write.csv(counts_pbmc, file.path(dir, "pbmc", "counts.csv"))
write.csv(coldata_pbmc, file.path(dir, "pbmc", "coldata.csv"))
write.csv(counts_tissue, file.path(dir, "tissue", "counts.csv"))
write.csv(coldata_tissue, file.path(dir, "tissue", "coldata.csv"))