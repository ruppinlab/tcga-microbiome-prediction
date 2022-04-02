suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggpubr)
  library(readr)
  library(VennDiagram)
})

color_palette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

top_hits <- read_tsv(
  "results/analysis/goodness_hits.txt",
  col_types = cols()
) %>%
  semi_join(
    read_tsv(
      "results/analysis/selected_hits.txt",
      col_types = cols()
    ),
    by = c("cancer", "analysis", "versus", "features", "how")
  ) %>%
  filter(analysis == "resp" & features %in% c("kraken", "htseq")) %>%
  select(cancer, versus, features, how, p_adj) %>%
  arrange(features, cancer, versus, p_adj)

# Take the hits from the top two methods
selected_hits <- top_hits %>%
  group_by(features, cancer, versus) %>%
  slice_head(n = 2) %>%
  mutate(rank = seq_len(n())) %>%
  ungroup() %>%
  select(cancer, versus, features, how, rank)

feature_table <- rbind(
  read_tsv(
    "results/analysis/microbial_features.txt",
    col_types = cols()
  ),
  read_tsv(
    "results/analysis/expression_features.txt",
    col_types = cols()
  )
)

# Filter the feature table to the selected hits.
feature_table <- feature_table %>%
  inner_join(selected_hits, by = c("cancer", "versus", "features", "how")) %>%
  select(cancer, versus, features, how, genera, median_rank, rank)

targets <- feature_table %>%
  select(cancer, versus, features) %>%
  unique()

firsts <- feature_table %>% filter(rank == 1)
seconds <- feature_table %>% filter(rank == 2)
mixed <- firsts %>%
  inner_join(seconds,
    by = c("cancer", "versus", "features", "genera"),
    suffix = c(".gold", ".silver")
  ) %>%
  select(
    cancer, versus, features, genera,
    how.gold, median_rank.gold, how.silver, median_rank.silver
  )

write_tsv(firsts, "figures/mlcomp/best_method.tsv")
write_tsv(seconds, "figures/mlcomp/second_best_method.tsv")
write_tsv(mixed, "figures/mlcomp/feature_correlation.tsv")

venn_summaries <- vector("list", nrow(targets))
for (i in seq_len(nrow(targets))) {
  cancer <- targets[i, "cancer", drop = TRUE]
  versus <- targets[i, "versus", drop = TRUE]
  features <- targets[i, "features", drop = TRUE]
  top_two_methods <- mixed %>%
    filter(cancer == !!cancer & versus == !!versus & features == !!features)

  gold <- top_two_methods[1, "how.gold", drop = TRUE]
  silver <- top_two_methods[1, "how.silver", drop = TRUE]
  filename <- paste0(
    tolower(cancer), "_", tolower(versus), "_", features, ".pdf"
  )

  correlation_plot <- ggscatter(
    top_two_methods,
    size = 4,
    x = "median_rank.gold",
    y = "median_rank.silver",
    add = "reg.line", conf.int = TRUE,
    cor.coef = TRUE, cor.method = "spearman",
    cor.coef.size = 8, # Affects the size of the cor coef
    cor.coeff.args = list(label.x.npc = "left", label.y.npc = "top"),
    xlab = paste(gold, "(median rank)"),
    ylab = paste(silver, "(median rank)"),
    title = paste(cancer, versus),
    font.title = 28,
    font.label = 28,
    font.x = 28,
    font.y = 28,
    add.params = list(color = "blue", fill = "lightgray"),
    font.tickslab = c(20, "plain", "black")
  )
  grid.newpage()
  pdf(paste0("figures/mlcomp/corr/", filename))
  print(correlation_plot)
  dev.off()

  grid.newpage()
  pdf(paste0("figures/mlcomp/venn/", filename))

  fill <- color_palette[c(3, 1)]
  if (features == "htseq") {
    fill <- color_palette[c(7, 1)]
  }

  venn_summary <- tibble(
    cancer, versus, features, gold, silver,
    area.gold = firsts %>%
      filter(cancer == !!cancer &
        versus == !!versus &
        features == !!features) %>%
      nrow(),
    area.silver = seconds %>%
      filter(cancer == !!cancer &
        versus == !!versus &
        features == !!features) %>%
      nrow(),
    co_occurance = top_two_methods %>% nrow()
  )

  venn_diagram <- draw.pairwise.venn(
    area1 = venn_summary$area.gold,
    area2 = venn_summary$area.silver,
    cross.area = venn_summary$co_occurance,
    c(venn_summary$gold, venn_summary$silver),
    fill = fill,
    fontface = "plain",
    fontfamily = "sans",
    cex = rep(3, 3),
    cat.cex = rep(3, 2),
    cat.fontfamily = "sans",
    cat.fontface = "plain",
    cat.default.pos = "outer",
    ind = FALSE
  )
  venn_figure <- cowplot::plot_grid(
    gTree(children = venn_diagram),
    scale = .93
  ) +
    cowplot::draw_label(
      label = paste(cancer, versus), x = 0.005, y = .995,
      hjust = 0, vjust = 1,
      fontfamily = "sans", fontface = "bold", size = 36
    )
  print(venn_figure)

  dev.off()
  venn_summaries[[i]] <- venn_summary
}
write_tsv(bind_rows(venn_summaries), "figures/mlcomp/venn_summary.tsv")
