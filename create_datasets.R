#!/usr/bin/env Rscript

options(warn=1)
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("GenomicDataCommons"))
suppressPackageStartupMessages(library("httr"))
suppressPackageStartupMessages(library("stringr"))

stopifnot(GenomicDataCommons::status()$status == "OK")

argp <- arg_parser("Create esets and other data")
argp <- add_argument(
    argp, "--data-dir", default="data", help="Data directory"
)
argp <- add_argument(
    argp, "--gdc-workflow-types", default=c("HTSeq - Counts"),
    help="GDC workflow types"
)
argp <- add_argument(
    argp, "--save-data-matrix", flag=TRUE,
    help="Save data matrix rds files"
)
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

if (args$save_data_matrix) {
    cat(
        "[", str_pad("Kraken Filtered", kraken_msg_pad, side="right"), "]",
        "Writing knight_kraken_meta.rds\n"
    )
    saveRDS(
        kraken_meta, paste(args$data_dir, "knight_kraken_meta.rds", sep="/")
    )
    cat(
        "[", str_pad("Kraken Filtered", kraken_msg_pad, side="right"), "]",
        "Writing knight_kraken_data.rds\n"
    )
    saveRDS(
        kraken_data, paste(args$data_dir, "knight_kraken_data.rds", sep="/")
    )
}

# GDC case metadata
days_per_year <- 365

gdc_cases_no_meta <- data.frame(
    case_uuid=c("375436b3-66ac-4d5e-b495-18a96d812a69"),
    case_submitter_id=c("TCGA-F5-6810")
)

cat("Getting GDC case metadata\n")
kraken_case_uuids <- unique(kraken_meta$case_uuid)
kraken_case_uuids <- kraken_case_uuids[
    !(kraken_case_uuids %in% gdc_cases_no_meta$case_uuid)
]
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
gdc_missing_case_meta$age_at_diagnosis <-
    gdc_missing_case_meta$age_at_diagnosis * days_per_year
gdc_case_meta <- rbind(gdc_case_meta, gdc_missing_case_meta)

gdc_case_meta$age_at_diagnosis[is.na(gdc_case_meta$age_at_diagnosis)] <-
    kraken_case_meta$age_at_diagnosis[
        kraken_case_meta$case_uuid
        %in% gdc_case_meta$case_uuid[is.na(gdc_case_meta$age_at_diagnosis)]
    ] * days_per_year
gdc_case_meta$tumor_stage <- gsub("^stage\\s+", "", gdc_case_meta$tumor_stage)
gdc_case_meta$tumor_stage <- gsub(
    "i/ii\\s+nos", "i or ii", gdc_case_meta$tumor_stage
)
gdc_case_meta$tumor_stage <- gsub("(a|b|c|s)$", "", gdc_case_meta$tumor_stage)
gdc_case_meta$tumor_stage[gdc_case_meta$tumor_stage == "not reported"] <- NA

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

# response metadata
cat("Processing drug response phenotypic data\n")
response_pdata <- read.delim(
    paste(args$data_dir, "drug_response.txt", sep="/"), row.names=NULL,
    stringsAsFactors=FALSE
)
colnames(response_pdata)[2] <- "case_submitter_id"
response_pdata$cancer <- paste("TCGA", response_pdata$cancers, sep="-")
response_pdata$cancers <- NULL
response_pdata <- response_pdata[
    c(ncol(response_pdata), 1:(ncol(response_pdata) - 1))
]
response_pdata[response_pdata == "[Not Available]"] <- NA
response_pdata[response_pdata == "[Unknown]"] <- NA
response_pdata$start.time[
    !is.na(response_pdata$start.time) &
    response_pdata$start.time == "[Completed]"
] <- 0
suppressWarnings(
    response_pdata$start.time <- as.integer(response_pdata$start.time)
)
response_pdata$end.time[
    !is.na(response_pdata$end.time) &
    response_pdata$end.time == "[Completed]"
] <- 0
suppressWarnings(
    response_pdata$end.time <- as.integer(response_pdata$end.time)
)
response_pdata <- response_pdata[order(
    response_pdata$cancer, response_pdata$case_submitter_id,
    response_pdata$drug.name, response_pdata$start.time
), ]
response_pdata <- response_pdata[
    !duplicated(response_pdata[c("case_submitter_id", "drug.name")]),
]
response_pdata <- merge(response_pdata, gdc_case_meta, by="case_submitter_id")
response_pdata$cancer <- as.factor(response_pdata$cancer)
response_pdata$response <- as.factor(response_pdata$response)
response_pdata$drug.name <- as.factor(response_pdata$drug.name)

# survival metadata
cat("Processing survival phenotypic data\n")
survival_pdata <- read.delim(
    paste(args$data_dir, "survival.txt", sep="/"), row.names=NULL,
    stringsAsFactors=FALSE
)
colnames(survival_pdata)[2] <- "case_submitter_id"
survival_pdata$cancer <- paste("TCGA", survival_pdata$cancer, sep="-")
survival_pdata <- survival_pdata[order(
    survival_pdata$cancer, survival_pdata$case_submitter_id
), ]
survival_pdata <- merge(survival_pdata, gdc_case_meta, by="case_submitter_id")

# generate datasets
min_uniq_cases <- 18
min_uniq_cases_per_class <- 4
min_uniq_case_exceptions <- c("stad oxaliplatin")

create_surv_eset <- function(
    adata, pdata, fdata, surv_type, data_type, msg_prefix,
    print_sample_msgs=TRUE
) {
    eset <- ExpressionSet(assayData=adata, phenoData=AnnotatedDataFrame(pdata))
    if (!is.null(fdata)) featureData(eset) <- AnnotatedDataFrame(fdata)
    status_col <- toupper(surv_type)
    time_col <- paste(status_col, "time", sep="_")
    eset$Status <- eset[[status_col]]
    eset$Survival_in_days <- eset[[time_col]]
    eset <- eset[
        , !is.na(eset[[status_col]]) & !is.na(eset[[time_col]])
        & eset[[time_col]] > 0
    ]
    pData(eset) <- gdata::drop.levels(pData(eset))
    if (anyDuplicated(eset$case_submitter_id)) {
        eset$Group <- match(
            eset$case_submitter_id, unique(eset$case_submitter_id)
        )
        if ("experimental_strategy" %in% colnames(pData(eset))) {
            case_strategy <- as.data.frame(table(data.frame(
                case_submitter_id=factor(eset$case_submitter_id),
                experimental_strategy=factor(eset$experimental_strategy)
            )))
            case_strategy <- case_strategy[case_strategy$Freq != 0, ]
            max_freq <- aggregate(Freq ~ case_submitter_id, case_strategy, max)
            colnames(max_freq)[colnames(max_freq) == "Freq"] <- "MaxFreq"
            case_strategy <-
                merge(case_strategy, max_freq, by="case_submitter_id")
            case_strategy$GroupWeight <-
                case_strategy$MaxFreq / case_strategy$Freq
            eset$GroupWeight <- case_strategy$GroupWeight[match(
                paste(eset$case_submitter_id, eset$experimental_strategy),
                paste(
                    case_strategy$case_submitter_id,
                    case_strategy$experimental_strategy
                )
            )]
        } else if ("GroupWeight" %in% colnames(pData(eset))) {
             eset$GroupWeight <- NULL
        }
    } else {
        if ("Group" %in% colnames(pData(eset))) eset$Group <- NULL
        if ("GroupWeight" %in% colnames(pData(eset))) eset$GroupWeight <- NULL
    }
    cancer_target <- gsub("^tcga ", "", paste(
        tolower(gsub("(-|\\s)+", " ", c(cancer, surv_type))), collapse=" "
    ))
    if (
        length(unique(eset$case_submitter_id)) >= min_uniq_cases
        || cancer_target %in% min_uniq_case_exceptions
    ) {
        if (print_sample_msgs) cat(
            msg_prefix, toupper(surv_type),
            length(unique(eset$case_submitter_id)), "cases",
            ncol(eset), "samples", "\n"
        )
        eset_filename_parts <- c(cancer, "surv", surv_type, data_type, "eset")
        eset_filename <- paste0(paste(
            tolower(gsub("(-|\\s)+", "_", eset_filename_parts)), collapse="_"
        ), ".rds")
        cat(msg_prefix, "Creating", eset_filename, "\n")
        saveRDS(eset, paste(args$data_dir, eset_filename, sep="/"))
    } else {
        cat(msg_prefix, paste("Skipping", toupper(surv_type)), "\n")
    }
}

create_drug_eset <- function(
    adata, pdata, fdata, resp_type, drug_name, data_type, msg_prefix,
    print_sample_msgs=TRUE
) {
    eset <- ExpressionSet(assayData=adata, phenoData=AnnotatedDataFrame(pdata))
    if (!is.null(fdata)) featureData(eset) <- AnnotatedDataFrame(fdata)
    eset$Class <- ifelse(
        eset$response %in% c(
            "Complete Response", "Partial Response"
        ), 1,
        ifelse(eset$response %in% c(
            "Clinical Progressive Disease", "Stable Disease"
        ), 0, NA)
    )
    eset$Class <- as.factor(eset$Class)
    eset <- eset[, !is.na(eset$Class)]
    pData(eset) <- gdata::drop.levels(pData(eset))
    if (anyDuplicated(eset$case_submitter_id)) {
        eset$Group <- match(
            eset$case_submitter_id, unique(eset$case_submitter_id)
        )
        num_uniq_cases_per_class <- aggregate(
            case_submitter_id ~ Class, pData(eset),
            function(x) length(unique(x))
        )$case_submitter_id
    } else {
        if ("Group" %in% colnames(pData(eset))) eset$Group <- NULL
        num_uniq_cases_per_class <- table(eset$Class)
    }
    cancer_target <- gsub("^tcga ", "", paste(
        tolower(gsub("(-|\\s)+", " ", c(cancer, drug_name))), collapse=" "
    ))
    if (
        (
            length(unique(eset$case_submitter_id)) >= min_uniq_cases
            || cancer_target %in% min_uniq_case_exceptions
        )
        && length(unique(eset$Class)) > 1
        && min(num_uniq_cases_per_class) >= min_uniq_cases_per_class
    ) {
        if (print_sample_msgs) cat(
            msg_prefix, drug_name, length(unique(eset$case_submitter_id)),
            "cases", ncol(eset), "samples\n"
        )
        eset_filename_parts <-
            c(cancer, resp_type, drug_name, data_type, "eset")
        eset_filename <- paste0(paste(
            tolower(gsub("(-|\\s)+", "_", eset_filename_parts)), collapse="_"
        ), ".rds")
        cat(msg_prefix, "Creating", eset_filename, "\n")
        saveRDS(eset, paste(args$data_dir, eset_filename, sep="/"))
    }
}

get_gdc_data <- function(project_id, workflow_types, msg_prefix) {
    file_query <-
        files() %>%
        filter(
            cases.project.project_id %in% project_id
            & analysis.workflow_type %in% workflow_types
        ) %>%
        GenomicDataCommons::select(c(
            "file_name",
            "analysis.workflow_type",
            "cases.project.project_id",
            "cases.case_id",
            "cases.submitter_id",
            "cases.samples.sample_id",
            "cases.samples.submitter_id",
            "cases.samples.sample_type",
            "cases.samples.is_ffpe",
            "cases.samples.portions.is_ffpe",
            "cases.samples.portions.analytes.aliquots.aliquot_id",
            "cases.samples.portions.analytes.aliquots.submitter_id"
        ))
    file_results <- results_all(file_query)
    file_meta <- data.frame(
        file_uuid=file_results$file_id,
        file_name=file_results$file_name,
        workflow_type=file_results$analysis$workflow_type,
        project_id=vapply(
            sapply(file_results$cases, `[[`, "project"), `[`, "project_id"
        ),
        case_uuid=sapply(file_results$cases, `[[`, "case_id"),
        case_submitter_id=sapply(file_results$cases, `[[`, "submitter_id"),
        sample_uuid=sapply(
            sapply(file_results$cases, `[[`, "samples"), `[[`, "sample_id"
        ),
        sample_submitter_id=sapply(
            sapply(file_results$cases, `[[`, "samples"), `[[`, "submitter_id"
        ),
        sample_type=sapply(
            sapply(file_results$cases, `[[`, "samples"), `[[`, "sample_type"
        ),
        sample_is_ffpe=sapply(
            sapply(file_results$cases, `[[`, "samples"), `[[`, "is_ffpe"
        ),
        portion_is_ffpe=sapply(
            sapply(
                sapply(file_results$cases, `[[`, "samples"), `[[`, "portions"
            ), `[[`, "is_ffpe"
        ),
        aliquot_uuid=sapply(
            sapply(
                sapply(
                    sapply(
                        sapply(
                            file_results$cases, `[[`, "samples"
                        ), `[[`, "portions"
                    ), `[[`, "analytes"
                ), `[[`, "aliquots"
            ), `[[`, "aliquot_id"
        ),
        aliquot_submitter_id=sapply(
            sapply(
                sapply(
                    sapply(
                        sapply(
                            file_results$cases, `[[`, "samples"
                        ), `[[`, "portions"
                    ), `[[`, "analytes"
                ), `[[`, "aliquots"
            ), `[[`, "submitter_id"
        ),
        row.names=file_results$file_id,
        stringsAsFactors=FALSE
    )
    file_meta <- file_meta[order(
        file_meta$project_id, file_meta$case_submitter_id,
        file_meta$sample_submitter_id, file_meta$aliquot_submitter_id
    ), ]
    cat(msg_prefix, "Downloading", length(file_results$id), "GDC files\n")
    files <- gdcdata(ids(file_query), progress=FALSE, use_cached=TRUE)
    return(list(meta=file_meta, files=files))
}

# GENCODE annots
gtf_annots_filename <- "gencode_v22_ensg_v98_annots.tsv"
gtf_annots_file <- paste(args$data_dir, gtf_annots_filename, sep="/")
cat("Loading", gtf_annots_file, "\n")
rna_annots <- read.delim(gtf_annots_file, row.names=1, stringsAsFactors=FALSE)
rna_annots$Symbol[grepl(".", rna_annots$Symbol, fixed=TRUE)] <- ""
combo_annots <- rbind(data.frame(
    Symbol=rep("", ncol(kraken_data)), row.names=colnames(kraken_data)
), rna_annots)

gdc_meta_pdata_cols <- c(
    "project_id", "case_uuid", "case_submitter_id", "sample_uuid",
    "sample_submitter_id", "sample_is_ffpe", "portion_is_ffpe", "aliquot_uuid",
    "aliquot_submitter_id", "file_uuid"
)
cancers <- sort(union(response_pdata$cancer, survival_pdata$cancer))
cancer_msg_pad <- max(str_length(cancers))
type_msg_pad <- max(str_length(c("kraken", "combo", args$gdc_workflow_types)))
surv_types <- c("os", "pfi")
resp_types <- c("resp")

cat("Generating datasets\n")
for (cancer in cancers) {
    # kraken
    data_type <- "kraken"
    msg_prefix <- paste("[",
        str_pad(cancer, cancer_msg_pad, side="right"),
        str_pad("Kraken", type_msg_pad, side="right"),
    "]")
    # survival
    kraken_surv_meta <- merge(
        survival_pdata[survival_pdata$cancer == cancer, , drop=FALSE],
        kraken_meta, by="case_submitter_id"
    )
    row.names(kraken_surv_meta) <- kraken_surv_meta$knight_id
    kraken_surv_meta$investigation <- NULL
    kraken_surv_meta <- kraken_surv_meta[c(2, 1, 3:ncol(kraken_surv_meta))]
    kraken_surv_data <-
        t(kraken_data)[, row.names(kraken_surv_meta), drop=FALSE]
    kraken_surv_meta <-
        kraken_surv_meta[order(row.names(kraken_surv_meta)), , drop=FALSE]
    kraken_surv_data <-
        kraken_surv_data[, order(colnames(kraken_surv_data)), drop=FALSE]
    for (surv_type in surv_types) {
        create_surv_eset(
            kraken_surv_data, kraken_surv_meta, NULL, surv_type, data_type,
            msg_prefix
        )
    }
    # response
    for (drug_name in sort(
        unique(response_pdata$drug.name[response_pdata$cancer == cancer])
    )) {
        kraken_drug_meta <- merge(
            response_pdata[
                response_pdata$cancer == cancer
                & response_pdata$drug.name == drug_name, , drop=FALSE
            ],
            kraken_meta, by="case_submitter_id"
        )
        if (nrow(kraken_drug_meta) == 0) next
        row.names(kraken_drug_meta) <- kraken_drug_meta$knight_id
        kraken_drug_meta$investigation <- NULL
        kraken_drug_meta <- kraken_drug_meta[c(2, 1, 3:ncol(kraken_drug_meta))]
        kraken_drug_data <-
            t(kraken_data)[, row.names(kraken_drug_meta), drop=FALSE]
        kraken_drug_meta <-
            kraken_drug_meta[order(row.names(kraken_drug_meta)), , drop=FALSE]
        kraken_drug_data <-
            kraken_drug_data[, order(colnames(kraken_drug_data)), drop=FALSE]
        for (idx in seq_along(resp_types)) {
            create_drug_eset(
                kraken_drug_data, kraken_drug_meta, NULL, resp_types[idx],
                drug_name, data_type, msg_prefix, idx == 1
            )
        }
    }
    # expression
    print_sample_msgs <- TRUE
    for (workflow_type in args$gdc_workflow_types) {
        msg_prefix <- paste("[",
            str_pad(cancer, cancer_msg_pad, side="right"),
            str_pad(workflow_type, type_msg_pad, side="right"),
        "]")
        gdc_data <- get_gdc_data(cancer, workflow_type, msg_prefix)
        # sample type filter
        gdc_data$meta <- gdc_data$meta[
            gdc_data$meta$sample_type %in% sample_types,
        ]
        cat(msg_prefix, "Generating data matrix\n")
        DT <- dcast(rbindlist(
            lapply(gdc_data$files[gdc_data$meta$file_uuid], fread),
            idcol="file_uuid"
        ), file_uuid ~ V1, value.var="V2")
        rna_data <-
            t(data.frame(DT[, !"file_uuid"], row.names=DT[["file_uuid"]]))
        rna_data <- rna_data[!grepl("^(X|N|_)_", row.names(rna_data)), ]
        rna_data <- rna_data[order(row.names(rna_data)), ]
        if (!identical(row.names(rna_annots), row.names(rna_data))) stop(
            msg_prefix, "RNA data matrix and annots row names not identical"
        )
        gdc_meta <- gdc_data$meta[, gdc_meta_pdata_cols]
        row.names(gdc_meta) <- gdc_meta$file_uuid
        if (args$save_data_matrix) {
            rna_data_filename_parts <- c(cancer, workflow_type)
            rna_data_filename <- paste0(paste(
                tolower(gsub("(-|\\s)+", "_", rna_data_filename_parts)),
                collapse="_"
            ), ".rds")
            cat(msg_prefix, "Writing", rna_data_filename, "\n")
            saveRDS(
                rna_data, paste(args$data_dir, rna_data_filename, sep="/")
            )
        }
        # survival
        rna_surv_meta <- merge(
            gdc_meta,
            survival_pdata[survival_pdata$cancer == cancer, , drop=FALSE],
            by="case_submitter_id"
        )
        row.names(rna_surv_meta) <- rna_surv_meta$file_uuid
        rna_surv_meta$cancer <- NULL
        rna_surv_meta <- rna_surv_meta[c(2, 3, 1, 4:ncol(rna_surv_meta))]
        rna_surv_data <- rna_data[, row.names(rna_surv_meta), drop=FALSE]
        rna_surv_meta <-
            rna_surv_meta[order(row.names(rna_surv_meta)), , drop=FALSE]
        rna_surv_data <-
            rna_surv_data[, order(colnames(rna_surv_data)), drop=FALSE]
        for (surv_type in surv_types) {
            create_surv_eset(
                rna_surv_data, rna_surv_meta, rna_annots, surv_type,
                workflow_type, msg_prefix, print_sample_msgs
            )
        }
        # combo survival
        combo_msg_prefix <- paste("[",
            str_pad(cancer, cancer_msg_pad, side="right"),
            str_pad("Combo", type_msg_pad, side="right"),
        "]")
        combo_surv_meta <- merge(
            rna_surv_meta,
            kraken_meta[!(colnames(kraken_meta) %in% c(
                "case_uuid", "case_submitter_id", "sample_submitter_id",
                "portion_is_ffpe"
            ))], by="sample_uuid"
        )
        row.names(combo_surv_meta) <- make.names(paste(
            combo_surv_meta$knight_id,
            ave(
                combo_surv_meta$knight_id, combo_surv_meta$knight_id,
                FUN=seq_along
            ), sep="_"
        ))
        combo_surv_meta <- combo_surv_meta[c(2:4, 1, 5:ncol(combo_surv_meta))]
        colnames(combo_surv_meta)[
            colnames(combo_surv_meta) == "aliquot_uuid.x"
        ] <- "gdc_aliquot_uuid"
        colnames(combo_surv_meta)[
            colnames(combo_surv_meta) == "aliquot_uuid.y"
        ] <- "knight_aliquot_uuid"
        colnames(combo_surv_meta)[
            colnames(combo_surv_meta) == "gdc_file_uuid"
        ] <- "knight_file_uuid"
        combo_surv_data <- rbind(
            t(kraken_data)[, combo_surv_meta$knight_id, drop=FALSE],
            rna_data[, combo_surv_meta$file_uuid, drop=FALSE]
        )
        colnames(combo_surv_data) <- row.names(combo_surv_meta)
        combo_surv_meta <-
            combo_surv_meta[order(row.names(combo_surv_meta)), , drop=FALSE]
        combo_surv_data <-
            combo_surv_data[, order(colnames(combo_surv_data)), drop=FALSE]
        if (!identical(row.names(combo_annots), row.names(combo_surv_data)))
            stop(
                msg_prefix,
                "Combo surv data matrix and annots row names not identical"
            )
        for (surv_type in surv_types) {
            create_surv_eset(
                combo_surv_data, combo_surv_meta, combo_annots, surv_type,
                "combo", combo_msg_prefix, print_sample_msgs
            )
        }
        # response
        for (drug_name in sort(
            unique(response_pdata$drug.name[response_pdata$cancer == cancer])
        )) {
            rna_drug_meta <- merge(
                gdc_meta,
                response_pdata[
                    response_pdata$cancer == cancer
                    & response_pdata$drug.name == drug_name, , drop=FALSE
                ], by="case_submitter_id"
            )
            if (nrow(rna_drug_meta) == 0) next
            row.names(rna_drug_meta) <- rna_drug_meta$file_uuid
            rna_drug_meta$cancer <- NULL
            rna_drug_meta <- rna_drug_meta[c(2, 3, 1, 4:ncol(rna_drug_meta))]
            rna_drug_data <- rna_data[, row.names(rna_drug_meta), drop=FALSE]
            rna_drug_meta <-
                rna_drug_meta[order(row.names(rna_drug_meta)), , drop=FALSE]
            rna_drug_data <-
                rna_drug_data[, order(colnames(rna_drug_data)), drop=FALSE]
            for (idx in seq_along(resp_types)) {
                create_drug_eset(
                    rna_drug_data, rna_drug_meta, rna_annots, resp_types[idx],
                    drug_name, workflow_type, msg_prefix,
                    ifelse(idx == 1, print_sample_msgs, FALSE)
                )
            }
            # combo response
            combo_drug_meta <- merge(
                rna_drug_meta,
                kraken_meta[!(colnames(kraken_meta) %in% c(
                    "case_uuid", "case_submitter_id", "sample_submitter_id",
                    "portion_is_ffpe"
                ))], by="sample_uuid"
            )
            row.names(combo_drug_meta) <- make.names(paste(
                combo_drug_meta$knight_id,
                ave(
                    combo_drug_meta$knight_id, combo_drug_meta$knight_id,
                    FUN=seq_along
                ), sep="_"
            ))
            combo_drug_meta <-
                combo_drug_meta[c(2:4, 1, 5:ncol(combo_drug_meta))]
            colnames(combo_drug_meta)[
                colnames(combo_drug_meta) == "aliquot_uuid.x"
            ] <- "gdc_aliquot_uuid"
            colnames(combo_drug_meta)[
                colnames(combo_drug_meta) == "aliquot_uuid.y"
            ] <- "knight_aliquot_uuid"
            colnames(combo_drug_meta)[
                colnames(combo_drug_meta) == "gdc_file_uuid"
            ] <- "knight_file_uuid"
            combo_drug_data <- rbind(
                t(kraken_data)[, combo_drug_meta$knight_id, drop=FALSE],
                rna_data[, combo_drug_meta$file_uuid, drop=FALSE]
            )
            colnames(combo_drug_data) <- row.names(combo_drug_meta)
            combo_drug_meta <-
                combo_drug_meta[order(row.names(combo_drug_meta)), , drop=FALSE]
            combo_drug_data <-
                combo_drug_data[, order(colnames(combo_drug_data)), drop=FALSE]
            if (!identical(row.names(combo_annots), row.names(combo_drug_data)))
                stop(
                    msg_prefix,
                    "Combo drug data matrix and annots row names not identical"
                )
            for (idx in seq_along(resp_types)) {
                create_drug_eset(
                    combo_drug_data, combo_drug_meta, combo_annots,
                    resp_types[idx], drug_name, "combo", combo_msg_prefix,
                    ifelse(idx == 1, print_sample_msgs, FALSE)
                )
            }
        }
        print_sample_msgs <- FALSE
    }
}
