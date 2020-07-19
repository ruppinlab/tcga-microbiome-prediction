suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
results_table <- read_tsv(args, col_types = cols())

# Filter out candidates that do not have good enough ROC/C-index.
# FDR correct the rest.
results_table <- results_table %>%
  filter(avg_test >= 0.6) %>%
  group_by(analysis, features) %>%
  mutate(p_adj = p.adjust(p_value, 'fdr')) %>%
  ungroup()

results_table %>%
  filter(p_greater <= 0.05 & p_adj <= 0.05) %>%
  arrange(analysis, features, desc(avg_test)) %>%
  mutate(
    avg_test = sprintf("%.3g", avg_test),
    sd_test = sprintf("%.3g", sd_test),
    avg_cov = sprintf("%.3g", avg_cov),
    p_adj = sprintf("%.2e", p_adj),
    p_value = sprintf("%.2e", p_value),
    p_greater = sprintf("%.2e", p_greater)
  ) %>%
  format_tsv() %>%
  cat()
