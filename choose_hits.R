suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
results_table <- read_tsv(args, col_types = cols())

results_table <- results_table %>%
  filter(!is.na(p_adj) & p_greater <= 0.05) %>%
  arrange(analysis, features, how, desc(avg_test)) %>%
  format_tsv() %>%
  cat()
