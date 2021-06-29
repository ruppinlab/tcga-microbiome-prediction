suppressMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
filename <- args[[1]]

features <- read_tsv(filename, col_types = cols())

stats <- features %>%
  group_by(cancer, what, features, how) %>%
  summarize(
    features = n(),
    pos = sum(mean > 0),
    neg = sum(mean < 0),
    min_seen = min(seen),
    first_q_seen = quantile(seen, probs = .25),
    median_seen = median(seen),
    third_q_seen = quantile(seen, probs = .75),
    max_seen = max(seen),
    .groups = "drop"
  )

stats %>%
  format_tsv() %>%
  cat()
