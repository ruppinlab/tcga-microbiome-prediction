options(warn=1)
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("httr"))
suppressPackageStartupMessages(library("GenomicDataCommons"))
suppressPackageStartupMessages(library("stringr"))

stopifnot(GenomicDataCommons::status()$status == "OK")

argp <- arg_parser("Process Knight kraken data and metadata")
argp <- add_argument(argp, "--data-dir", default="data", help="Data directory")
argp <- add_argument(argp, "--debug", flag=TRUE, help="Debug mode")
args <- parse_args(argp)

# Knight Kraken data
kraken_msg_pad <- 19
sample_msg_pad <- 5

cat("Processing Knight Kraken data\n")
kraken_base_url <-
    "ftp://ftp.microbio.me/pub/cancer_microbiome_analysis/TCGA/Kraken/"
kraken_data_filename <- "Kraken-TCGA-Voom-SNM-Plate-Center-Filtering-Data.csv"
kraken_data_url <- paste0(kraken_base_url, kraken_data_filename)
kraken_data_file <- paste(args$data_dir, kraken_data_filename, sep="/")
if (file.exists(kraken_data_file)) {
    cat("Using existing", kraken_data_file, "\n")
} else {
    cat("Downloading", kraken_data_filename, "\n")
    download.file(kraken_data_url, kraken_data_file)
}
kraken_data <- as.matrix(fread(kraken_data_file), rownames=1)
kraken_meta_filename <- "Metadata-TCGA-Kraken-17625-Samples.csv"
kraken_meta_url <- paste0(kraken_base_url, kraken_meta_filename)
kraken_meta_file <- paste(args$data_dir, kraken_meta_filename, sep="/")
if (file.exists(kraken_meta_file)) {
    cat("Using existing", kraken_meta_file, "\n")
} else {
    cat("Downloading", kraken_meta_filename, "\n")
    download.file(kraken_meta_url, kraken_meta_file)
}
kraken_meta <- read.csv(
    kraken_meta_file, row.names=NULL, stringsAsFactors=FALSE
)
colnames(kraken_meta)[1] <- "knight_id"
row.names(kraken_meta) <- kraken_meta$knight_id
kraken_meta$case_uuid <- tolower(kraken_meta$case_uuid)
kraken_meta$sample_uuid <- tolower(kraken_meta$sample_uuid)
kraken_meta$aliquot_uuid <- tolower(kraken_meta$aliquot_uuid)
kraken_meta$gdc_file_uuid <- tolower(kraken_meta$gdc_file_uuid)
kraken_meta <- kraken_meta[order(
    kraken_meta$investigation, kraken_meta$experimental_strategy,
    kraken_meta$case_uuid, kraken_meta$sample_uuid,
    kraken_meta$aliquot_uuid, kraken_meta$data_submitting_center_label
), ]
kraken_covar_cols <- c("age_at_diagnosis", "gender", "pathologic_stage_label")
cat(
    "[", str_pad("Kraken", kraken_msg_pad, side="right"), "]",
    nrow(kraken_meta), "samples", length(unique(kraken_meta$case_uuid)),
    "unique cases\n"
)

if (args$debug) {
    kraken_meta <- kraken_meta[, c(
        "investigation", "case_uuid", kraken_covar_cols, "knight_id",
        "sample_uuid", "aliquot_uuid", "filename", "experimental_strategy",
        "sample_type", "reference_genome", "data_submitting_center_label",
        "portion_is_ffpe"
    )]
}

# sample type filter
sample_types <- c("Primary Tumor", "Additional - New Primary")
kraken_meta <- kraken_meta[kraken_meta$sample_type %in% sample_types, ]

# WGS filters
wgs_main_centers <- c(
    "Baylor College of Medicine",
    "Broad Institute of MIT and Harvard",
    "Washington University School of Medicine"
)
wgs_meta <- kraken_meta[kraken_meta$experimental_strategy == "WGS", ]
cat(
    "[", str_pad("Kraken WGS", kraken_msg_pad, side="right"), "]",
    str_pad(nrow(wgs_meta), sample_msg_pad, side="left"), "samples",
    str_pad(length(unique(wgs_meta$case_uuid)), sample_msg_pad, side="left"),
    "unique cases\n"
)
wgs_case_freq <- table(factor(wgs_meta$case_uuid))
wgs_dupes <- wgs_meta[
    wgs_meta$case_uuid %in% names(wgs_case_freq)[wgs_case_freq > 1],
]
wgs_dupes_no_main_only_case_uuids <- setdiff(
    wgs_dupes$case_uuid[
        !(wgs_dupes$data_submitting_center_label %in% wgs_main_centers)
    ],
    wgs_dupes$case_uuid[
        wgs_dupes$data_submitting_center_label %in% wgs_main_centers
    ]
)
wgs_meta <- wgs_meta[
    wgs_meta$case_uuid %in% names(wgs_case_freq)[wgs_case_freq == 1]
    | wgs_meta$case_uuid %in% wgs_dupes_no_main_only_case_uuids
    | wgs_meta$data_submitting_center_label
      == "Broad Institute of MIT and Harvard"
    | (
        wgs_meta$data_submitting_center_label
        == "Washington University School of Medicine"
        & grepl("GRCh37-lite", wgs_meta$reference_genome)
    )
    | (
        wgs_meta$data_submitting_center_label == "Baylor College of Medicine"
        & (
            wgs_meta$reference_genome == "GRCh37-lite"
            | grepl("Illumina", wgs_meta$filename, ignore.case=TRUE)
        )
    ),
]
cat(
    "[", str_pad("Kraken WGS Filtered", kraken_msg_pad, side="right"), "]",
    str_pad(nrow(wgs_meta), sample_msg_pad, side="left"), "samples",
    str_pad(length(unique(wgs_meta$case_uuid)), sample_msg_pad, side="left"),
    "unique cases\n"
)

# RNA-seq filters
rna_meta <- kraken_meta[kraken_meta$experimental_strategy == "RNA-Seq", ]
cat(
    "[", str_pad("Kraken RNA", kraken_msg_pad, side="right"), "]",
    str_pad(nrow(rna_meta), sample_msg_pad, side="left"), "samples",
    str_pad(length(unique(rna_meta$case_uuid)), sample_msg_pad, side="left"),
    "unique cases\n"
)
rna_case_freq <- table(factor(rna_meta$case_uuid))
rna_meta <- rna_meta[
    rna_meta$case_uuid %in% names(rna_case_freq)[rna_case_freq == 1]
    | (
        rna_meta$data_submitting_center_label
        == "Broad Institute of MIT and Harvard"
        & rna_meta$reference_genome == "HG19_Broad_variant"
    )
    | (
        rna_meta$data_submitting_center_label
        == "University of North Carolina"
        & grepl("sorted_genome_alignments", rna_meta$filename)
    )
    | (
        rna_meta$data_submitting_center_label
        == "Canada's Michael Smith Genome Sciences Centre"
        & rna_meta$reference_genome == "GRCh37-lite"
    ),
]
cat(
    "[", str_pad("Kraken RNA Filtered", kraken_msg_pad, side="right"), "]",
    str_pad(nrow(rna_meta), sample_msg_pad, side="left"), "samples",
    str_pad(length(unique(rna_meta$case_uuid)), sample_msg_pad, side="left"),
    "unique cases\n"
)

kraken_meta <- rbind(wgs_meta, rna_meta)
kraken_meta <- kraken_meta[order(row.names(kraken_meta)), ]
kraken_data <- kraken_data[row.names(kraken_meta), ]
cat(
    "[", str_pad("Kraken Filtered", kraken_msg_pad, side="right"), "]",
    str_pad(nrow(kraken_meta), sample_msg_pad, side="left"), "samples",
    str_pad(length(unique(kraken_meta$case_uuid)), sample_msg_pad, side="left"),
    "unique cases\n"
)

# GDC case metadata
days_per_year <- 365.2422
gdc_cases_no_meta <- data.frame(
    case_uuid=c("375436b3-66ac-4d5e-b495-18a96d812a69"),
    case_submitter_id=c("TCGA-F5-6810")
)

kraken_case_uuids <- unique(kraken_meta$case_uuid)
kraken_case_uuids <- kraken_case_uuids[
    !(kraken_case_uuids %in% gdc_cases_no_meta$case_uuid)
]
gdc_case_results_filename <- "gdc_case_results.rds"
gdc_case_results_file <- paste(
    args$data_dir, gdc_case_results_filename, sep="/"
)
if (file.exists(gdc_case_results_file)) {
    cat("Using existing", gdc_case_results_file, "\n")
    gdc_case_results <- readRDS(gdc_case_results_file)
} else {
    cat("Downloading GDC case metadata\n")
    gdc_case_query <-
        cases() %>%
        filter(case_id %in% kraken_case_uuids) %>%
        GenomicDataCommons::select(c(
            "submitter_id",
            "diagnoses.age_at_diagnosis",
            "demographic.gender",
            "diagnoses.tumor_stage"
        ))
    gdc_case_results <- results_all(gdc_case_query)
    cat(paste("Writing", gdc_case_results_filename), "\n")
    saveRDS(gdc_case_results, gdc_case_results_file)
}
# set NULL diagnoses inner list elements to NA
for (i in seq_along(gdc_case_results$diagnoses)) {
    if (is.null(gdc_case_results$diagnoses[[i]]))
        gdc_case_results$diagnoses[[i]] <- data.frame(
            age_at_diagnosis=NA, tumor_stage=NA
        )
}
gdc_case_meta <- data.frame(
    case_uuid=gdc_case_results$case_id,
    case_submitter_id=gdc_case_results$submitter_id,
    gender=gdc_case_results$demographic$gender,
    age_at_diagnosis=sapply(
        gdc_case_results$diagnoses, `[[`, "age_at_diagnosis"
    ),
    tumor_stage=sapply(gdc_case_results$diagnoses, `[[`, "tumor_stage"),
    row.names=NULL,
    stringsAsFactors=FALSE
)
kraken_case_meta <- kraken_meta[, c("case_uuid", kraken_covar_cols)]
colnames(kraken_case_meta)[
    colnames(kraken_case_meta) == "pathologic_stage_label"
] <- "tumor_stage"
kraken_case_meta <- kraken_case_meta[!duplicated(kraken_case_meta), ]
gdc_missing_case_meta <- merge(
    gdc_cases_no_meta, kraken_case_meta, by="case_uuid"
)
gdc_missing_case_meta$gender <- tolower(gdc_missing_case_meta$gender)
gdc_missing_case_meta$tumor_stage <- tolower(gdc_missing_case_meta$tumor_stage)
gdc_missing_case_meta$age_at_diagnosis <- as.integer(round(
    gdc_missing_case_meta$age_at_diagnosis * days_per_year
))
gdc_case_meta <- rbind(gdc_case_meta, gdc_missing_case_meta)

gdc_case_meta$age_at_diagnosis[is.na(gdc_case_meta$age_at_diagnosis)] <-
    as.integer(round(kraken_case_meta$age_at_diagnosis[
        match(gdc_case_meta$case_uuid, kraken_case_meta$case_uuid)
    ] * days_per_year))[is.na(gdc_case_meta$age_at_diagnosis)]
gdc_case_meta$tumor_stage <- gsub("^stage\\s+", "", gdc_case_meta$tumor_stage)
gdc_case_meta$tumor_stage <-
    gsub("i/ii\\s+nos", "i or ii", gdc_case_meta$tumor_stage)
gdc_case_meta$tumor_stage <- gsub("(a|b|c|s)$", "", gdc_case_meta$tumor_stage)
gdc_case_meta$tumor_stage[
    gdc_case_meta$tumor_stage %in% c("not reported", "x")
] <- NA

kraken_meta$case_submitter_id <- gdc_case_meta$case_submitter_id[
    match(kraken_meta$case_uuid, gdc_case_meta$case_uuid)
]
kraken_meta <- kraken_meta[, c(
    "case_uuid", "case_submitter_id", "knight_id",
    colnames(kraken_meta)[!(
        colnames(kraken_meta)
        %in% c("case_uuid", "case_submitter_id", "knight_id")
    )]
)]
gdc_case_meta$case_uuid <- NULL

kraken_meta <- kraken_meta[, !(colnames(kraken_meta) %in% kraken_covar_cols)]

gdc_aliquot_results_filename <- "gdc_aliquot_results.rds"
gdc_aliquot_results_file <- paste(
    args$data_dir, gdc_aliquot_results_filename, sep="/"
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
    aliquot_uuid=unlist(
        gdc_aliquot_results$aliquot_ids, use.names=FALSE
    ),
    aliquot_submitter_id=unlist(
        gdc_aliquot_results$submitter_aliquot_ids, use.names=FALSE
    ),
    row.names=NULL,
    stringsAsFactors=FALSE
)
kraken_meta$aliquot_submitter_id <- gdc_aliquot_meta$aliquot_submitter_id[
    match(kraken_meta$aliquot_uuid, gdc_aliquot_meta$aliquot_uuid)
]

cat("Writing knight_kraken_meta.rds\n")
saveRDS(kraken_meta, paste(args$data_dir, "knight_kraken_meta.rds", sep="/"))
cat("Writing knight_kraken_data.rds\n")
saveRDS(kraken_data, paste(args$data_dir, "knight_kraken_data.rds", sep="/"))
cat("Writing gdc_case_meta.rds\n")
saveRDS(gdc_case_meta, paste(args$data_dir, "gdc_case_meta.rds", sep="/"))
