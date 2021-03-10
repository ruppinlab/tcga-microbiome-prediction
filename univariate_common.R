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


read_kraken_aliquots <- function(metafile, mapfile) {
  meta <- readRDS("data/knight_kraken_meta.rds")
  aliquot_map <- read_tsv("data/aliquot_map.tsv", col_types = cols())

  tibble(
    sample_id = rownames(meta), aliquot_uuid = meta$aliquot_uuid
  ) %>%
    inner_join(aliquot_map, by = "aliquot_uuid") %>%
    select(sample_id, aliquot_barcode) %>%
    mutate(sample_barcode = aliquot_to_sample(aliquot_barcode))
}

find_microbial_sample <- function(kraken, aliquots, j) {
  kraken[, j] %>%
    enframe(name = "sample_id", value = "abundance") %>%
    inner_join(aliquots, by = "sample_id") %>%
    select(aliquot_barcode, sample_barcode, abundance)
}


join_response_with_microbial <- function(these_res, microbial_sample) {
  d <- these_res %>%
    select(sample_barcode = case_submitter_id, response) %>%
    inner_join(microbial_sample, by = "sample_barcode")
  d <- d[sample(seq_len(nrow(d))), ]
  d <- d %>%
    filter(!duplicated(sample_barcode)) %>%
    arrange(sample_barcode)
  d
}


join_survival_with_microbial <- function(surv, microbial_sample) {
  data <- surv %>%
    select(sample_barcode = case_submitter_id, survival) %>%
    inner_join(microbial_sample, by = "sample_barcode")
  data <- data[sample(seq_len(nrow(data))), ]
  data <- data %>% filter(!duplicated(sample_barcode))
  data
}
