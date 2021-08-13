library(dplyr)
library(readr)
library("ggpubr")
library(VennDiagram)
library(cowplot)

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

features <- rbind(
  read_tsv("results/analysis/microbial_features.txt"),
  read_tsv("results/analysis/expression_features.txt")
)

selected_hits <- read_tsv("results/analysis/goodness_hits.txt") %>%
  semi_join(read_tsv("results/analysis/selected_hits.txt")) %>%
  filter(analysis == "resp" & features %in% c("kraken", "htseq")) %>%
  select(cancer, versus, features, how, avg_test) %>%
  arrange(features, cancer, versus, desc(avg_test)) %>%
  group_by(features, cancer, versus) %>%
  slice_head(n = 2) %>%
  mutate(rank = seq(1, 2)) %>%
  ungroup() %>%
  select(cancer, versus, features, how, rank)

feats <- features %>%
  inner_join(selected_hits, by = c("cancer", "versus", "features", "how")) %>%
  select(cancer, versus, features, how, genera, median_rank, rank)

targets <- feats %>%
  select(cancer, versus, features) %>%
  unique()

firsts <- feats %>% filter(rank == 1)
write_tsv(firsts, "figures/mlcomp/best_method.tsv")
seconds <- feats %>% filter(rank == 2)
write_tsv(firsts, "figures/mlcomp/second_best_method.tsv")
mixed <- firsts %>%
  inner_join(seconds, by = c("cancer", "versus", "features", "genera"))
write_tsv(firsts, "figures/mlcomp/joined_methods.tsv")

for (i in seq_len(nrow(targets))) {
  cancer <- targets[i, "cancer", drop = TRUE]
  versus <- targets[i, "versus", drop = TRUE]
  features <- targets[i, "features", drop = TRUE]
  my_data <- mixed %>%
    filter(cancer == !!cancer & versus == !!versus & features == !!features)

  area.x <- firsts %>%
    filter(cancer == !!cancer & versus == !!versus & features == !!features) %>%
    nrow()
  area.y <- seconds %>%
    filter(cancer == !!cancer & versus == !!versus & features == !!features) %>%
    nrow()
  cross <- my_data %>% nrow()
  filename <- paste0(
    tolower(cancer), "_", tolower(versus), "_", features, ".pdf"
  )

  pl <- ggscatter(my_data,
    x = "median_rank.x", y = "median_rank.y",
    add = "reg.line", conf.int = TRUE,
    cor.coef = TRUE, cor.method = "spearman",
    cor.coef.size = 5.5, # Affects the size of the cor coef
    xlab = paste(my_data[1, "how.x", drop = TRUE], "(median rank)"),
    ylab = paste(my_data[1, "how.y", drop = TRUE], "(median rank)"),
    title = paste(cancer, versus),
    font.title = 18,
    font.label = 18,
    font.x = 16,
    font.y = 16,
    add.params = list(color = "blue", fill = "lightgray")
  )
  grid.newpage()
  pdf(paste0("figures/mlcomp/corr/", filename))
  print(pl)
  dev.off()

  grid.newpage()
  pdf(paste0("figures/mlcomp/venn/", filename))

  fill <- cbPalette[c(3, 1)]
  if (features == "htseq") {
    fill <- cbPalette[c(7, 1)]
  }

  g <- draw.pairwise.venn(
    area.x,
    area.y,
    cross,
    c(my_data[1, "how.x", drop = TRUE], my_data[1, "how.y", drop = TRUE]),
    fill = fill,
    fontface = "bold",
    fontfamily = "sans",
    cex = rep(1.5, 3),
    cat.cex = rep(1.5, 2),
    cat.fontfamily = "sans",
    cat.fontface = "bold",
    cat.default.pos = "outer",
    ind = FALSE
  )
  p <- cowplot::plot_grid(gTree(children = g)) +
    draw_label(
      label = paste(cancer, versus), x = .2, y = .95,
      fontfamily = "sans", fontface = "bold", size = 20
    )
  print(p)

  dev.off()
}
