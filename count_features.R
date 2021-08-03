suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

features <- read_tsv("results/analysis/simplified_microbial.txt",
  col_types = cols()
)
colnames(features) <- sub(" ", "_", colnames(features))
features <- features %>%
  mutate(Univariate_FDR = as.numeric(
    ifelse(Univariate_FDR == "n.s.", "1", Univariate_FDR)
  ))

cat("Genera appeared\n")
genera_counts <- features %>%
  group_by(Genus) %>%
  tally() %>%
  select(Genus, n)
cat("Once", sum(genera_counts$n == 1))
cat("\n")
cat("Twice", sum(genera_counts$n == 2))
cat("\n")
cat("More than twice", sum(genera_counts$n > 2))
cat("\n")
cat("Total", nrow(genera_counts))
cat("\n\n")

cat("####\n")
features %>%
  filter(Univariate_FDR <= 0.05) %>%
  pull(Genus) %>%
  unique() %>%
  length() %>%
  cat()
cat(" features are univariately predictive in some model")

cat("\n\n")

cat("####\n")

features %>%
  group_by(Cancer, Versus) %>%
  tally() %>%
  format_tsv() %>%
  cat()
cat("\n\n")

cat("####\n")

cat("Distribution of features per model\n")
features %>%
  group_by(Cancer, Versus) %>%
  tally() %>%
  pull(n) %>%
  summary() %>%
  print()
cat("\n\n")

cat("####\n")

kingdoms <- features %>%
  distinct(Genus) %>%
  mutate(kingdom = sapply(strsplit(Genus, "[.]"), function(x) x[[1]])) %>%
  group_by(kingdom) %>%
  tally()
kingdoms <- rbind(kingdoms, tibble(kingdom = "All", n = sum(kingdoms$n)))
kingdoms %>%
  format_tsv() %>%
  cat()
cat("\n\n")

cat("####\n")

cat("Phyla with multiplicity\n")
phylum <- features %>%
  select(Genus) %>%
  mutate(
    phylum = sapply(
      strsplit(Genus, "[.]"),
      function(x) paste(x[c(1, 2)], collapse = ".")
    )
  ) %>%
  group_by(phylum) %>%
  tally(name = "phylum_with_multiplicity") %>%
  arrange(desc(phylum_with_multiplicity))

phylum %>%
  format_tsv() %>%
  cat()
cat("\n\n")

cat("####\n")
inconsistent_directions <-
  features %>%
  select(Conditional_Direction) %>%
  filter(!(Conditional_Direction %in% c(
    "Positive", "Positive,Positive", "Positive,Positive,Positive",
    "Negative", "Negative,Negative", "Negative,Negative,Negative"
  )))
cat(nrow(inconsistent_directions), "had inconsistent directions\n")

features <- features %>%
  mutate(
    Conditional_Direction = ifelse(
      grepl("^Positive", Conditional_Direction),
      "Positive",
      "Negative"
    )
  )

cat("\n\n")

cat("####\n")
all_microbes <- features %>%
  group_by(Conditional_Direction) %>%
  tally(name = "microbe_count")

all_bacteria <- features %>%
  filter(startsWith(Genus, "k__Bacteria")) %>%
  group_by(Conditional_Direction) %>%
  tally(name = "bacteria_count")

firmicutes <- features %>%
  filter(startsWith(Genus, "k__Bacteria.p__Firmicutes")) %>%
  group_by(Conditional_Direction) %>%
  tally(name = "firmicutes_count")

joint_count <- all_microbes %>%
  inner_join(all_bacteria, by = "Conditional_Direction") %>%
  inner_join(firmicutes, by = "Conditional_Direction")

joint_count %>%
  format_tsv() %>%
  cat()
cat("\n\n####\n")

print(
  fisher.test(
    as.matrix(joint_count[, c("bacteria_count", "firmicutes_count")])
  )
)

cat("\n\n####\n")

binom.test(
  joint_count[1, "microbe_count", drop = TRUE],
  sum(joint_count[1:2, "microbe_count"])
)

cat("\n\n####\n")

cat("Firmicute counts")
features %>%
  filter(startsWith(Genus, "k__Bacteria.p__Firmicutes")) %>%
  group_by(Cancer, Conditional_Direction) %>%
  tally(name = "Count") %>%
  pivot_wider(
    names_from = "Conditional_Direction", values_from = "Count",
    values_fill = 0
  ) %>%
  mutate(Total = Negative + Positive) %>%
  rename(
    `Selected Negatively` = Negative, `Selected Positively` =
      Positive
  ) %>%
  format_tsv() %>%
  cat()
