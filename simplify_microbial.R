suppressPackageStartupMessages({
  library(tidyverse)
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
      univariate_direction = ifelse(
        univariate_direction == -1, "negative", "positive"
      ),
      univariate_fdr = sprintf("%.3g", univariate_fdr)
    ) %>%
    select(
      cancer,
      versus,
      genus = genera,
      median_rank,
      models_present = seen,
      direction,
      univariate_direction,
      univariate_fdr
    )
}

features <- bind_rows(results)

features %>%
  arrange(cancer, versus, median_rank, desc(models_present)) %>%
  format_tsv() %>%
  cat()
