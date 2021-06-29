library(tidyverse)

vs_features <- list()

for (file in list.files("response_vs_features",
  pattern = "[.]tsv$",
  full.names = TRUE
)) {
  vs_features[[length(vs_features) + 1]] <- read_tsv(file, col_types = cols())
}

for (file in list.files("survival_vs_features",
  pattern = "[.]tsv$",
  full.names = TRUE
)) {
  vs_features[[length(vs_features) + 1]] <- read_tsv(file, col_types = cols())
}

vs_features <- bind_rows(vs_features)
vs_features <- vs_features %>%
  group_by(cancer, what) %>%
  mutate(p_bh = p.adjust(p_value, "BH")) %>%
  ungroup()

feature_hits <- vs_features %>%
  filter(p_bh <= 0.05) %>%
  arrange(cancer, what, p_bh)
cat(format_tsv(feature_hits))
