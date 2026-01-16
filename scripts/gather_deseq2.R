# Capture output and error messages to log file
log <- file(snakemake@log[[1]], open="wt")
sink(log, type="output")
sink(log, type="message")

library(readr)
library(purrr)

res_combined <- snakemake@input %>%
    map_df(~read_csv(.x))

write_csv(res_combined, snakemake@output[[1]])