suppressMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
features <- args[1]
hits <- read_tsv(args[2], col_types = cols())
hits <- hits %>%
  filter(features == !!features) %>%
  mutate(features = ifelse(
    features == "kraken", features, paste0(features, "_counts")
  ))

cancer <- tolower(hits$cancer)
versus <- tolower(hits$versus)
how <- tolower(hits$how)
name <- paste(
  "tcga", cancer, hits$analysis, versus, hits$features, how,
  sep = "_"
)
path <- file.path(
  "results/models",
  hits$analysis,
  name,
  paste(name, "feature_weights.rds", sep = "_")
)
cat(path, sep = "\n")
