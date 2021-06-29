suppressMessages({
  library(dplyr)
  library(readr)
})

args = commandArgs(trailingOnly = TRUE)
hits <- read_tsv(args[1], col_types = cols())
hits <- hits %>%
  filter(features == "kraken" &
     ( ( analysis == "resp" & how == "RFE" ) |
       ( analysis == "surv" & how == "CNET" ) ) )
 
cancer <- tolower(hits$cancer)
versus <- tolower(hits$versus)
how <- tolower(hits$how)
name <- paste(
  "tcga", cancer, hits$analysis, versus, hits$features, how, sep = "_"
)
path <- file.path(
  "results/models",
  hits$analysis,
  name,
  paste(name, "feature_weights.rds", sep = "_")
)
cat(path, sep = "\n")
