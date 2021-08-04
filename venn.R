suppressPackageStartupMessages({
  library(VennDiagram)
  library(RColorBrewer)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
microbial_features_txt <- args[[1]]
expression_features_txt <- args[[2]]
figures_venn <- args[[3]]

myCol <- brewer.pal(3, "Pastel2")

microbial_features <- read_tsv(microbial_features_txt, col_types = cols()) %>%
  filter(what == "resp") %>%
  mutate(path = paste(cancer, versus, genera, sep = "."))
microbial_sets <- split(
  microbial_features$path,
  factor(microbial_features$how, levels = c("LGR", "LIMMA", "RFE"))
)

venn.diagram(microbial_sets,
  filename = file.path(figures_venn, "microbial_features.png"),
  imagetype = "png",
  fill = myCol,
  cex = 1.0,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.2,
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.default.pos = "outer"
)

expression_features <- read_tsv(expression_features_txt, col_types = cols()) %>%
  filter(what == "resp") %>%
  mutate(path = paste(cancer, versus, genera, sep = "."))
expression_sets <- split(
  expression_features$path,
  factor(expression_features$how, levels = c("LGR", "EDGER", "RFE"))
)

grid.newpage()
venn.diagram(expression_sets,
  filename = file.path(figures_venn, "expression_features.png"),
  imagetype = "png",
  fill = myCol,
  cex = 1.0,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.2,
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.default.pos = "outer"
)
