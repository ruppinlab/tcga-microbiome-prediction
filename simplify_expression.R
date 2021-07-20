suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

results <- vector("list", length(args))
for (k in seq_along(args)) {
  results[[k]] <- read_tsv(args[[k]], col_types = cols()) %>%
    filter((what == "surv" & how == "CNET") |
      (what == "resp" & how == "RFE")) %>%
    mutate(
      direction = ifelse(
        p_value > 0.05, NA, ifelse(p_greater > .1, "negative", "positive")
      ),
      symbol = sapply(
        mapIds(org.Hs.eg.db,
          keys = sub("[.].*$", "", genera),
          column = "SYMBOL", keytype = "ENSEMBL", multiVals = "list"
        ),
        paste,
        collapse = ","
      )
    ) %>%
    select(
      cancer,
      versus,
      symbol,
      ensembl = genera,
      median_rank,
      models_present = seen,
      direction
    ) %>%
    mutate(symbol = ifelse(symbol == "NA", NA, symbol))
}

features <- bind_rows(results)

features %>%
  arrange(cancer, versus, median_rank, desc(models_present)) %>%
  format_tsv() %>%
  cat()
