suppressMessages({
  library(dplyr)
  library(readr)
})

hits <- read_tsv("analysis/goodness_hits.txt", col_types = cols())
hits <- hits %>% filter(features == "kraken" & analysis %in% c("resp", "surv"))

cancer <- tolower(hits$cancer)
versus <- tolower(hits$versus)
how <- ifelse(hits$analysis == "surv", "cnet", "rfe")
name <- paste(
  "tcga", cancer, hits$analysis, versus, hits$features, how, sep = "_"
)
path <- file.path(
  "results",
  hits$analysis,
  name,
  paste(name, "feature_weights.rds", sep = "_")
)
cat(path, sep = "\n")
