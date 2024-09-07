options(warn = 1)
suppressPackageStartupMessages({
    library(argparser)
    library(httr)
    library(GenomicDataCommons)
    library(missForest)
    library(stringr)
})

set.seed(777)

stopifnot(GenomicDataCommons::status()$status == "OK")

argp <- arg_parser("Process new WGS Kraken2+Bracken data and metadata")
argp <- add_argument(
    argp, "--data-dir",
    default = "data", help = "Data directory"
)
argp <- add_argument(argp, "--debug", flag = TRUE, help = "Debug mode")
args <- parse_args(argp)

cat("Processing new WGS Kraken2+Bracken data\n")

kraken_data_filename <- "tcga_wgs_primary_tumors_genus_count_matrix.tsv"
kraken_data_file <- paste(args$data_dir, kraken_data_filename, sep = "/")
cat("Loading", kraken_data_filename, "\n")
kraken_data <- read.delim(
    kraken_data_file,
    sep = "\t", header = TRUE, check.names = FALSE
)
kraken_data <- kraken_data[
    !(kraken_data$taxonomy_id %in% c(9605, 9606)),
]
kraken_feature_meta <- kraken_data[c("taxonomy_id", "name", "taxonomy_lvl")]
row.names(kraken_feature_meta) <- kraken_feature_meta$name
row.names(kraken_data) <- kraken_data$name
kraken_feature_meta$name <- NULL
kraken_data$taxonomy_id <- NULL
kraken_data$name <- NULL
kraken_data$taxonomy_lvl <- NULL
kraken_data <- as.matrix(kraken_data)
storage.mode(kraken_data) <- "integer"

kraken_sample_meta_filename <- "tcga_wgs_primary_tumors_file_meta.tsv"
kraken_sample_meta_file <- paste(
    args$data_dir, kraken_sample_meta_filename,
    sep = "/"
)
cat("Loading", kraken_sample_meta_filename, "\n")
kraken_sample_meta <- read.delim(
    kraken_sample_meta_file,
    sep = "\t", header = TRUE, check.names = FALSE
)
row.names(kraken_sample_meta) <- kraken_sample_meta$file_id

# GDC case metadata
gdc_case_results_filename <- "gdc_case_results.rds"
gdc_case_results_file <- paste(
    args$data_dir, gdc_case_results_filename,
    sep = "/"
)
if (file.exists(gdc_case_results_file)) {
    cat("Using existing", gdc_case_results_file, "\n")
    gdc_case_results <- readRDS(gdc_case_results_file)
} else {
    cat("Downloading GDC case metadata\n")
    gdc_case_query <-
        cases() %>%
        GenomicDataCommons::filter(
            case_id %in% kraken_sample_meta$case_id
        ) %>%
        GenomicDataCommons::select(c(
            "project.project_id",
            "submitter_id",
            "demographic.gender",
            "diagnoses.age_at_diagnosis",
            "diagnoses.ajcc_pathologic_stage"
        ))
    gdc_case_results <- results_all(gdc_case_query)
    cat(paste("Writing", gdc_case_results_filename), "\n")
    saveRDS(gdc_case_results, gdc_case_results_file)
}
for (i in seq_along(gdc_case_results$diagnoses)) {
    # set NULL diagnoses inner list elements to NA
    if (is.null(gdc_case_results$diagnoses[[i]])) {
        gdc_case_results$diagnoses[[i]] <- data.frame(
            age_at_diagnosis = NA, ajcc_pathologic_stage = NA
        )
    } else if (
        # set non-existent diagnoses.ajcc_pathologic_stage to "Not Reported"
        is.null(gdc_case_results$diagnoses[[i]]$ajcc_pathologic_stage)
    ) {
        gdc_case_results$diagnoses[[i]]$ajcc_pathologic_stage <- "Not Reported"
    }
}
gdc_case_meta <- data.frame(
    project_id = gdc_case_results$project$project_id,
    case_id = gdc_case_results$case_id,
    case_submitter_id = gdc_case_results$submitter_id,
    gender = gdc_case_results$demographic$gender,
    age_at_diagnosis = sapply(
        gdc_case_results$diagnoses, `[[`, "age_at_diagnosis"
    ),
    tumor_stage = sapply(
        gdc_case_results$diagnoses, `[[`, "ajcc_pathologic_stage"
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
)
gdc_case_meta$gender <- tolower(gdc_case_meta$gender)
gdc_case_meta$gender[
    gdc_case_meta$project_id %in% c("TCGA-BRCA", "TCGA-OV")
        & is.na(gdc_case_meta$gender)
] <- "female"
gdc_case_meta$gender[
    gdc_case_meta$project_id %in% c("TCGA-TGCT")
        & is.na(gdc_case_meta$gender)
] <- "male"
gdc_case_meta$gender[is.na(gdc_case_meta$gender)] <- "unknown"
gdc_case_meta$tumor_stage <- tolower(gdc_case_meta$tumor_stage)
gdc_case_meta$tumor_stage <- gsub(
    "^stage\\s+", "", gdc_case_meta$tumor_stage
)
gdc_case_meta$tumor_stage <- gsub(
    "i/ii\\s+nos", "i or ii", gdc_case_meta$tumor_stage
)
gdc_case_meta$tumor_stage <- gsub(
    "(a|b|c|s)$", "", gdc_case_meta$tumor_stage
)
gdc_case_meta$tumor_stage <- gsub(
    "^(not reported|x)$", NA, gdc_case_meta$tumor_stage
)

xmis <- gdc_case_meta[
    c("project_id", "gender", "age_at_diagnosis", "tumor_stage")
]
xmis$project_id <- factor(xmis$project_id)
xmis$gender <- factor(xmis$gender)
xmis$tumor_stage <- ordered(
    xmis$tumor_stage,
    levels = c("i", "ii", "iii", "iv")
)
imp <- missForest(xmis)
gdc_case_meta$age_at_diagnosis <- as.integer(round(imp$ximp$age_at_diagnosis))
gdc_case_meta$tumor_stage <- as.character(imp$ximp$tumor_stage)

cat("Writing k2b_kraken_sample_meta.rds\n")
saveRDS(
    kraken_sample_meta,
    paste(args$data_dir, "k2b_kraken_sample_meta.rds", sep = "/")
)
cat("Writing k2b_kraken_feature_meta.rds\n")
saveRDS(
    kraken_feature_meta,
    paste(args$data_dir, "k2b_kraken_feature_meta.rds", sep = "/")
)
cat("Writing k2b_kraken_data.rds\n")
saveRDS(
    kraken_data, paste(args$data_dir, "k2b_kraken_data.rds", sep = "/")
)
cat("Writing gdc_case_meta.rds\n")
saveRDS(gdc_case_meta, paste(args$data_dir, "gdc_case_meta.rds", sep = "/"))
