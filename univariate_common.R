suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

aliquot_to_sample <- function(aliquot_barcode) {
  sapply(
    strsplit(aliquot_barcode, "-"),
    function(x) paste0(x[1:3], collapse = "-")
  )
}


read_kraken <- function(file) {
  readRDS(file)
}


read_kraken_aliquots <- function(metafile) {
  meta <- readRDS(metafile)
  as_tibble(meta, rownames = "sample_id") %>%
    select(sample_id,
      aliquot_barcode = aliquot_submitter_id,
      sample_barcode = case_submitter_id
    )
}

find_microbial_sample <- function(kraken, aliquots, j) {
  kraken[, j] %>%
    enframe(name = "sample_id", value = "abundance") %>%
    inner_join(aliquots, by = "sample_id") %>%
    select(aliquot_barcode, sample_barcode, abundance)
}


join_response_with_microbial <- function(these_res, microbial_sample) {
  d <- these_res %>%
    inner_join(microbial_sample, by = "sample_barcode") %>%
    select(sample_barcode, abundance, response)
  # Randomly permute the rows so that a unique, but random, row is chosen
  d <- d[sample(seq_len(nrow(d))), ]
  d %>%
    filter(!duplicated(sample_barcode)) %>%
    arrange(sample_barcode)
}


join_survival_with_microbial <- function(surv, microbial_sample) {
  data <- surv %>%
    select(sample_barcode = case_submitter_id, survival) %>%
    inner_join(microbial_sample, by = "sample_barcode")
  data <- data[sample(seq_len(nrow(data))), ]
  data <- data %>% filter(!duplicated(sample_barcode))
  data
}
