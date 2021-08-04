suppressPackageStartupMessages({
  library(tidyverse)
})

collapse_conditional_direction <- function(direction) {
  if (length(direction) == 1) {
    return(direction)
  } else if (all(direction == direction[[1]])) {
    return(paste0(direction[[1]], '*'))
  } else {
    return(paste(direction, collapse=','))
  }
}

args <- commandArgs(trailingOnly = TRUE)
feature_files <- args[1:2]
microbial_features <- read_tsv(args[[3]], col_types = cols())

results <- vector("list", length(feature_files))
for (k in seq_along(feature_files)) {
  hits <-
    read_tsv(feature_files[[k]], col_types = cols()) %>%
    inner_join(
      microbial_features,
      by = c("cancer", "what", "versus", "features", "genera")
    )
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
      univariate_fdr = ifelse(is.na(univariate_fdr) | univariate_fdr > 0.05,
        "n.s.",
        sprintf("%.3g", univariate_fdr)
      )
    ) %>%
    select(
      Cancer = cancer,
      Versus = versus,
      Genus = genera,
      How = how,
      `Models Present` = seen,
      `Conditional Direction` = direction,
      `Univariate FDR` = univariate_fdr,
      `Univariate Direction` = univariate_direction
    )
}

features <- bind_rows(results)
features <- features %>%
  group_by(Cancer, Versus, Genus) %>%
  summarize(
    How = paste(How, collapse = ","),
    `Models Present` = paste0(sum(`Models Present`), '/', n() * 100),
    `Conditional Direction` = collapse_conditional_direction(`Conditional Direction`),
    `Univariate FDR` = `Univariate FDR`[1],
    `Univariate Direction` = `Univariate Direction`[1],
    .groups = "drop"
  )

features %>%
  arrange(Cancer, Versus, desc(`Models Present`)) %>%
  format_tsv() %>%
  cat()
