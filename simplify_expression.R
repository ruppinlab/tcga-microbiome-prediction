suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
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

selected_hits <- read_tsv(args[[2]], col_types=cols()) %>%
  select(cancer, what=analysis, versus, features, how)

features <- read_tsv(args[[1]], col_types = cols()) %>%
  semi_join(selected_hits, by=c('cancer', 'what', 'versus', 'features', 'how'))


hits <- features %>% group_by(cancer, what, versus, features, genera)%>%
  tally(name = 'method_count') %>%
  ungroup %>%
  filter(what == 'surv' | method_count >= 2)

features <- features %>%
  semi_join(hits, by=c('cancer', 'what', 'versus', 'features', 'genera')) %>%
  mutate(conditional_direction = ifelse(p_greater <= 0.5, 'Positive', 'Negative')) %>%
  select(cancer, versus, how, genera, seen, conditional_direction)

features <- features %>%
  group_by(cancer, versus, genera) %>%
  summarize(
    seen = paste0(sum(seen), '/', n() * 100),
    how = paste(how, collapse=','),
    conditional_direction = collapse_conditional_direction(conditional_direction),
    .groups = 'drop')

features <- features %>% 
  mutate(
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
    How = how,
    `Models Present` = seen,
    `Conditional Direction` = conditional_direction
  ) %>%
  mutate(`Gene Symbol` = ifelse(`Gene Symbol` == "NA", NA, `Gene Symbol`))

features %>%
  arrange(Cancer, Versus, desc(`Models Present`)) %>%
  format_tsv() %>%
  cat()
