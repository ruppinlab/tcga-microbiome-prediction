options(warn = 1)
suppressPackageStartupMessages({
    library(argparser)
    library(httr)
    library(GenomicDataCommons)
    library(readxl)
    library(stringr)
})

stopifnot(GenomicDataCommons::status()$status == "OK")

argp <- arg_parser("Process Salzberg kraken data and metadata")
argp <- add_argument(
    argp, "--data-dir",
    default = "data", help = "Data directory"
)
argp <- add_argument(argp, "--debug", flag = TRUE, help = "Debug mode")
args <- parse_args(argp)

# Salzberg Kraken data
kraken_msg_pad <- 15
sample_msg_pad <- 5

cat("Processing Salzberg Kraken data\n")

# kraken_base_url <-
#     "https://github.com/yge15/TCGA_Microbial_Content/raw/main/"
# kraken_data_filename <-
#     "TableS1_Kraken-TCGA-WGS-5734-Samples-11349-Species.Microbial2023_noEuk.RawCounts.xlsx"
# kraken_data_url <- paste0(kraken_base_url, kraken_data_filename)
# kraken_data_file <- paste(args$data_dir, kraken_data_filename, sep = "/")
# if (file.exists(kraken_data_file)) {
#     cat("Using existing", kraken_data_file, "\n")
# } else {
#     cat("Downloading", kraken_data_filename, "\n")
#     download.file(kraken_data_url, kraken_data_file)
# }
# kraken_data <- as.data.frame(
#     read_excel(kraken_data_file, progress = FALSE),
#     stringsAsFactors = FALSE
# )
# row.names(kraken_data) <- kraken_data$HopkinsID
# kraken_data$HopkinsID <- NULL
# kraken_data <- as.matrix(kraken_data)
# storage.mode(kraken_data) <- "integer"
# kraken_data <- kraken_data[!(rownames(kraken_data) == "Sum"), ]

kraken_data <- readRDS(
    paste(args$data_dir, "salzberg_kraken_data.rds", sep = "/")
)
kraken_data <- as.matrix(kraken_data)
storage.mode(kraken_data) <- "integer"

# kraken_meta_filename <- "TableS12_Metadata-TCGA-WGS-5734-Samples.xlsx"
# kraken_meta_url <- paste0(kraken_base_url, kraken_meta_filename)
# kraken_meta_file <- paste(args$data_dir, kraken_meta_filename, sep = "/")
# if (file.exists(kraken_meta_file)) {
#     cat("Using existing", kraken_meta_file, "\n")
# } else {
#     cat("Downloading", kraken_meta_filename, "\n")
#     download.file(kraken_meta_url, kraken_meta_file)
# }
# kraken_meta <- as.data.frame(
#     read_excel(kraken_meta_file, progress = FALSE),
#     stringsAsFactors = FALSE
# )
# colnames(kraken_meta) <- c(
#     "hopkins_id", "knight_id", "staussman_id", "file_uuid",
#     "project_id", "case_uuid", "sample_uuid", "sample_type",
#     "aliquot_uuid", "aliquot_submitter_id", "rd_len", "tot_rds",
#     "grch38_unmapped", "grch38+chm13_unmapped"
# )
# row.names(kraken_meta) <- kraken_meta$hopkins_id
# kraken_meta$case_uuid <- tolower(kraken_meta$case_uuid)
# kraken_meta$sample_uuid <- tolower(kraken_meta$sample_uuid)
# kraken_meta$aliquot_uuid <- tolower(kraken_meta$aliquot_uuid)
# kraken_meta$file_uuid <- tolower(kraken_meta$file_uuid)

kraken_meta <- readRDS(
    paste(args$data_dir, "salzberg_kraken_meta.rds", sep = "/")
)

cat(
    "[", str_pad("Kraken", kraken_msg_pad, side = "right"), "]",
    nrow(kraken_meta), "samples", length(unique(kraken_meta$case_uuid)),
    "unique cases\n"
)

# # sample type filter
# sample_types <- c(
#     "primary_tumor", "primary_blood_derived_cancer_-_peripheral_blood"
# )
# kraken_meta <- kraken_meta[kraken_meta$sample_type %in% sample_types, ]

# kraken_data <- kraken_data[row.names(kraken_meta), ]
# cat(
#     "[", str_pad("Kraken Filtered", kraken_msg_pad, side = "right"), "]",
#     nrow(kraken_meta), "samples", length(unique(kraken_meta$case_uuid)),
#     "unique cases\n"
# )

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
        filter(case_id %in% kraken_meta$case_uuid) %>%
        GenomicDataCommons::select(c(
            "project.project_id",
            "submitter_id",
            "diagnoses.age_at_diagnosis",
            "demographic.gender",
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
    case_uuid = gdc_case_results$case_id,
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
kraken_meta$case_submitter_id <- gdc_case_meta$case_submitter_id[
    match(kraken_meta$case_uuid, gdc_case_meta$case_uuid)
]
gdc_case_meta$case_uuid <- NULL

gdc_aliquot_results_filename <- "gdc_aliquot_results.rds"
gdc_aliquot_results_file <- paste(
    args$data_dir, gdc_aliquot_results_filename,
    sep = "/"
)
if (file.exists(gdc_aliquot_results_file)) {
    cat("Using existing", gdc_aliquot_results_file, "\n")
    gdc_aliquot_results <- readRDS(gdc_aliquot_results_file)
} else {
    cat("Downloading GDC aliquot metadata\n")
    gdc_aliquot_query <-
        cases() %>%
        filter(
            samples.portions.analytes.aliquots.aliquot_id
            %in% kraken_meta$aliquot_uuid
        ) %>%
        GenomicDataCommons::select(c(
            "aliquot_ids",
            "submitter_aliquot_ids"
        ))
    gdc_aliquot_results <- results_all(gdc_aliquot_query)
    cat(paste("Writing", gdc_aliquot_results_filename), "\n")
    saveRDS(gdc_aliquot_results, gdc_aliquot_results_file)
}
gdc_aliquot_meta <- data.frame(
    aliquot_uuid = unlist(
        gdc_aliquot_results$aliquot_ids,
        use.names = FALSE
    ),
    aliquot_submitter_id = unlist(
        gdc_aliquot_results$submitter_aliquot_ids,
        use.names = FALSE
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
)
kraken_meta$aliquot_submitter_id <- gdc_aliquot_meta$aliquot_submitter_id[
    match(kraken_meta$aliquot_uuid, gdc_aliquot_meta$aliquot_uuid)
]

cat("Writing salzberg_kraken_meta.rds\n")
saveRDS(
    kraken_meta, paste(args$data_dir, "salzberg_kraken_meta.rds", sep = "/")
)
cat("Writing salzberg_kraken_data.rds\n")
saveRDS(
    kraken_data, paste(args$data_dir, "salzberg_kraken_data.rds", sep = "/")
)
cat("Writing gdc_case_meta.rds\n")
saveRDS(gdc_case_meta, paste(args$data_dir, "gdc_case_meta.rds", sep = "/"))
cat("Writing gdc_aliquot_meta.rds\n")
saveRDS(gdc_aliquot_meta, paste(args$data_dir, "gdc_aliquot_meta.rds", sep = "/"))
