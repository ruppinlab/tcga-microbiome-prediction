suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

do_wilcox <- function(x, y) {
  test <- NULL
  try({
    test <- wilcox.test(x, y, paired = TRUE)
  })
  ifelse(is.null(test), as.double(NA), test$p.value)
}

do_onesided_wilcox <- function(x, y) {
  test <- NULL
  try({
    test <- wilcox.test(x, y, paired = TRUE, alt = "g")
  })
  ifelse(is.null(test), as.double(NA), test$p.value)
}

args <- commandArgs(trailingOnly = TRUE)
tests <- read_tsv(args[1], col_types = cols())
covariates <- read_tsv(args[2], col_types = cols())

covariates <- bind_rows(list(
  covariates,
  # Limma and Edger compare to LGR.  Duplicate the covariates *but*
  # remember that LIMMA is only kraken and Edger only expression.
  # The unmatched rows will be removed in the inner join.
  covariates %>%
    filter(how == "LGR") %>%
    mutate(how = "LIMMA"),
  covariates %>% filter(how == "LGR") %>% mutate(how = "EDGER")
))

joined_goodness <- inner_join(
  covariates,
  tests,
  by = c("cancer", "analysis", "versus", "features", "how", "index"),
  suffix = c("_cov", "_test"),
)

results_table <-
  joined_goodness %>%
  group_by(cancer, analysis, versus, features, how) %>%
  summarize(
    avg_test = mean(goodness_test, na.rm = TRUE),
    sd_test = sd(goodness_test, na.rm = TRUE),
    avg_cov = mean(goodness_cov, na.rm = TRUE),
    sd_cov = sd(goodness_cov, na.rm = TRUE),
    p_value = ifelse(
      avg_test >= 0.6, do_wilcox(goodness_test, goodness_cov), NA
    ),
    p_greater = ifelse(
      avg_test >= 0.6, do_onesided_wilcox(goodness_test, goodness_cov), NA
    ),
    .groups = "drop"
  )

adjusted <- results_table %>%
  filter(!is.na(p_value)) %>%
  group_by(analysis, features, how) %>%
  mutate(p_adj = p.adjust(p_value, "fdr")) %>%
  ungroup()

not_adjusted <- results_table %>%
  filter(is.na(p_value)) %>%
  mutate(p_adj = NA)

rbind(adjusted, not_adjusted) %>%
  arrange(analysis, features, how, desc(avg_test)) %>%
  format_tsv() %>%
  cat()
