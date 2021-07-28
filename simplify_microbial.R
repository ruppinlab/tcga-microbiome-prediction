suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
feature_files <- args[1:2]
selected_hits <- read_tsv(args[[3]], col_types=cols()) %>%
  select(cancer, what=analysis, versus, features, how)

results <- vector("list", length(feature_files))
for (k in seq_along(feature_files)) {
  hits <- read_tsv(feature_files[[k]], col_types = cols()) %>%
    semi_join(selected_hits, by=c('cancer', 'what', 'versus', 'features', 'how'))

   results[[k]] <- hits %>%
     mutate(
      direction = ifelse(
        p_value > 0.05, NA, ifelse(p_greater > .1, "Negative", "Positive")
      ),
      univariate_direction =
        ifelse(
          is.na(univariate_fdr) | univariate_fdr > 0.05,
          "",
          ifelse(
            univariate_direction == -1, "Negative", "Positive"
          )
        ),
      univariate_fdr = ifelse(is.na(univariate_fdr) | univariate_fdr > 0.05, 'n.s.', sprintf("%.3g", univariate_fdr))
     ) %>%
     select(
      Cancer = cancer,
      Versus = versus,
      Genus = genera,
      `Median Rank` = median_rank,
      `Models Present` = seen,
      `Conditional Direction` = direction,
      `Univariate FDR` = univariate_fdr,
      `Univariate Direction` = univariate_direction
    )
}

features <- bind_rows(results)

features %>%
  arrange(Cancer, Versus, `Median Rank`, desc(`Models Present`)) %>%
  format_tsv() %>%
  cat()
