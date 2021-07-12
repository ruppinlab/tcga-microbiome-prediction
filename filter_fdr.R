suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
cutoff <- as.numeric(args[1])
results_table <- read_tsv(args[2], col_types = cols())

results_table %>%
  filter(p_adj <= cutoff) %>%
  arrange(analysis, features, how, desc(avg_test)) %>%
  mutate(
    avg_test = sprintf("%.4g", avg_test),
    sd_test = sprintf("%.4g", sd_test),
    avg_cov = sprintf("%.4g", avg_cov),
    sd_cov = sprintf("%.4g", sd_cov),
    p_adj = sprintf("%.2e", p_adj),
    p_value = sprintf("%.2e", p_value),
    p_greater = sprintf("%.2e", p_greater)
  ) %>%
  format_tsv() %>%
  cat()
