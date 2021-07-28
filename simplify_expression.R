suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

selected_hits <- read_tsv(args[[2]], col_types=cols()) %>%
  select(cancer, what=analysis, versus, features, how)

features <- read_tsv(args[[1]], col_types = cols()) %>%
  semi_join(selected_hits, by=c('cancer', 'what', 'versus', 'features', 'how'))

features <- features %>% 
  mutate(
    direction = ifelse(
      p_value > 0.05, NA, ifelse(p_greater > .1, "Negative", "Positive")
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
    Cancer = cancer,
    Versus = versus,
    `Gene Symbol` = symbol,
    ENSEMBL = genera,
    `Median Rank`=median_rank,
    `Models Present` = seen,
    `Conditional Direction` = direction
  ) %>%
  mutate(`Gene Symbol` = ifelse(`Gene Symbol` == "NA", NA, `Gene Symbol`))

features %>%
  arrange(Cancer, Versus, `Median Rank`, desc(`Models Present`)) %>%
  format_tsv() %>%
  cat()
