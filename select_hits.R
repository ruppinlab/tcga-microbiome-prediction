# Select CNET survival models unconditionally and RFE response models
# if at least one other method validates the model
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
hits <- read_tsv(args[[1]], col_types=cols())
ml_counts <- hits %>%
  group_by(cancer, analysis, versus, features) %>%
  tally(name = "ml_counts")

hits <- hits %>%
  inner_join(
    ml_counts,
    by = c("cancer", "analysis", "versus", "features")
  ) %>%
  filter(
    (analysis == "surv" & how == "CNET") |
      (analysis == "resp" &  ml_counts > 1)
  ) %>%
  select(cancer, analysis, versus, features, how) %>%
  arrange(features, analysis, cancer, versus)

hits %>%
  format_tsv() %>%
  cat()
