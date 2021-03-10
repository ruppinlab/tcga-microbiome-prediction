suppressPackageStartupMessages({
  library(tidyverse)
  library(coin)
})
source("univariate_common.R")

args <- commandArgs(trailingOnly = TRUE)
cancer <- args[[1]]
long_cancer <- paste0('TCGA-', cancer)
what <- args[[2]]
what_time <- paste0(what, "_time")

features <- read_tsv("analysis/microbial_features.txt", col_types = cols())

set.seed(98765)
# A seed for each iteration
eff_seeds <- sample(1:2^15, 10000)
k <- 1

surv <- readRDS("data/survival_pdata.rds") %>%
  select(case_submitter_id, cancer, s = one_of(what), s_time = one_of(what_time)) %>%
  filter(cancer == long_cancer & !is.na(s_time))

mid <- median(Surv(surv$s_time, surv$s), na.rm = TRUE)
mid <- as.numeric(mid[["quantile"]])

# Remove persons lost to follow up before median time
surv <- surv %>% filter(s_time >= mid | s == 1)
surv <- surv %>% mutate(
  survival = as.factor(ifelse(s_time >= mid, "b_long", "a_short"))
)

kraken <- read_kraken("data/knight_kraken_data.rds")
aliquots <- read_kraken_aliquots(
  "data/knight_kraken_meta.rds",
  "data/aliquot_map.tsv"
)
results <- list()
genera <- features %>%
  filter(cancer == !!cancer, what == !!what) %>%
  pull(genera) %>%
  unique()
 
for (genus in genera) {
  j <- which(colnames(kraken) == genus)
  microbial_sample <- find_microbial_sample(kraken, aliquots, j) %>%
    select(sample_barcode, abundance)

  seed <- eff_seeds[k]
  k <- k + 1

  set.seed(seed)
  data <- join_survival_with_microbial(surv, microbial_sample)

  if (length(unique(data$survival)) != 2) {
    next
  }

  test <- NULL
  p_value <- NULL
  direction <- 0
  try({
    wc <- wilcox_test(
      formula = abundance ~ survival, data = data,
      distribution = "approximate"
    )
    if (TRUE || pvalue(wc) < 0.01) {
      wc <- wilcox_test(
        formula = abundance ~ survival, data = data,
        distribution = "exact"
      )
    }
    if (pvalue(wc) < 0.01) {
      ranks <- rank(data$abundance)
      mean_rank_affected <- mean(ranks[data$survival == "b_long"])
      mean_rank_unaffected <- mean(ranks[data$survival == "a_short"])

      if (mean_rank_affected > mean_rank_unaffected) {
        direction <- 1
      } else {
        direction <- -1
      }
    }
    p_value <- pvalue(wc)
  })
  if (!is.numeric(p_value)) {
    next
  }

  results[[length(results) + 1]] <-
    tibble(
      cancer = cancer, what = what, genus = genus,
      direction = direction, p_value = p_value, seed = seed
    )
}

results <- do.call(rbind, results) %>% arrange(p_value)

cat(format_tsv(results))
