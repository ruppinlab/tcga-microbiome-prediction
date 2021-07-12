suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
results_table <- read_tsv(args, col_types = cols())

# Filter out candidates that do not have good enough ROC/C-index.
# FDR correct the rest.
adjusted <- results_table %>%
  filter(!is.na(p_value)) %>%
  group_by(analysis, features, how) %>%
  mutate(p_adj = p.adjust(p_value, "fdr")) %>%
  ungroup()

not_adjusted <- results_table %>%
  filter(is.na(p_value)) %>%
  mutate(p_adj = NA)

rbind(adjusted, not_adjusted) %>%
  arrange(analysis, features, how, desc(avg_test)) %>%
  format_tsv() %>%
  cat()
