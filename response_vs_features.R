suppressPackageStartupMessages({
  library(tidyverse)
  library(coin)
})
source("do_cibersort_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
cancer <- args[[1]]
drug_name <- args[[2]]

features <- read_tsv("analysis/microbial_features.txt", col_types = cols())

set.seed(98765)
# A seed for each iteration
eff_seeds <- sample(1:2^15, 10000)
k <- 1

res <- readRDS("data/response_pdata.rds") %>%
  select(case_submitter_id, cancer,
    drug_name = drug.name, start_time = start.time,
    response
  ) %>%
  mutate(
    response = as.factor(
      ifelse(response %in% c("Complete Response", "Partial Response"),
        "no", "yes"
      )
    ),
    cancer = sub("TCGA-", "", cancer)
  )

kraken <- read_kraken("data/knight_kraken_data.rds")
aliquots <- read_kraken_aliquots(
  "data/knight_kraken_meta.rds",
  "data/aliquot_map.tsv"
)
results <- list()
genera <- features %>%
  filter(cancer == !!cancer, what == !!drug_name) %>%
  pull(genera) %>%
  unique()

these_res <- res %>%
  filter(cancer == !!cancer & drug_name == !!drug_name) %>%
  arrange(case_submitter_id, cancer, start_time)

for (genus in genera) {
  j <- which(colnames(kraken) == genus)
  microbial_sample <- find_microbial_sample(kraken, aliquots, j) %>%
    select(sample_barcode, abundance)

  seed <- eff_seeds[k]
  k <- k + 1

  set.seed(seed)
  data <- join_response_with_microbial(these_res, microbial_sample)

  if (length(unique(data$response)) != 2) {
    next
  }

  test <- NULL
  p_value <- NULL
  direction <- 0
  try({
    wc <- wilcox_test(
      formula = abundance ~ response, data = data,
      distribution = "approximate"
    )
    if (TRUE || pvalue(wc) < 0.01) {
      wc <- wilcox_test(
        formula = abundance ~ response, data = data,
        distribution = "exact"
      )
    }
    if (pvalue(wc) < 0.01) {
      ranks <- rank(data$abundance)
      mean_rank_affected <- mean(ranks[data$response == "yes"])
      mean_rank_unaffected <- mean(ranks[data$response == "no"])

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
      cancer = cancer, what = drug_name, genus = genus,
      direction = direction, p_value = p_value, seed = seed
    )
}

results <- do.call(rbind, results) %>% arrange(p_value)

cat(format_tsv(results))
