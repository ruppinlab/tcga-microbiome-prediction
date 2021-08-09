suppressPackageStartupMessages(library("Biobase"))

uniq_kraken_case_uuids <- c()
uniq_kraken_na_case_uuids <- c()

for (eset_file in Sys.glob("data/tcga_*_kraken_eset.rds")) {
    cat("Loading", eset_file, "\n")
    eset <- readRDS(eset_file)
    pdata <- pData(eset)
    uniq_kraken_case_uuids <- union(uniq_kraken_case_uuids, pdata$case_uuid)
    uniq_kraken_na_case_uuids <- union(
        uniq_kraken_na_case_uuids, pdata$case_uuid[is.na(pdata$tumor_stage)]
    )
}

cat("Kraken cases:", length(uniq_kraken_case_uuids), "\n")
cat(
    "Kraken NA tumor stage cases:", length(uniq_kraken_na_case_uuids),
    "\n\n"
)

uniq_rna_case_uuids <- c()
uniq_rna_na_case_uuids <- c()

for (eset_file in Sys.glob("data/tcga_*_htseq_counts_eset.rds")) {
    cat("Loading", eset_file, "\n")
    eset <- readRDS(eset_file)
    pdata <- pData(eset)
    uniq_rna_case_uuids <- union(uniq_rna_case_uuids, pdata$case_uuid)
    uniq_rna_na_case_uuids <- union(
        uniq_rna_na_case_uuids, pdata$case_uuid[is.na(pdata$tumor_stage)]
    )
}

cat("RNA cases:", length(uniq_rna_case_uuids), "\n")
cat("RNA NA tumor stage cases:", length(uniq_rna_na_case_uuids), "\n")
