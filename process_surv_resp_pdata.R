options(warn=1)
suppressPackageStartupMessages(library(argparser))

argp <- arg_parser("Process survival and drug response phenotypic data")
argp <- add_argument(
    argp, "--data-dir", default="data", help="Data directory"
)
args <- parse_args(argp)

cat("Loading gdc_case_meta.rds\n")
gdc_case_meta <- readRDS(paste(args$data_dir, "gdc_case_meta.rds", sep="/"))

# response metadata
cat("Processing drug response phenotypic data\n")
response_pdata <- read.delim(
    paste(args$data_dir, "tcga_drug_response.tsv", sep="/"), row.names=NULL,
    stringsAsFactors=FALSE
)
response_pdata$cancer <- paste("TCGA", response_pdata$cancer, sep="-")
response_pdata[response_pdata == "[Not Available]"] <- NA
response_pdata[response_pdata == "[Discrepancy]"] <- NA
response_pdata[response_pdata == "[Unknown]"] <- NA
response_pdata$start_time[
    !is.na(response_pdata$start_time) &
    response_pdata$start_time == "[Completed]"
] <- 0
suppressWarnings(
    response_pdata$start_time <- as.integer(response_pdata$start_time)
)
response_pdata$end_time[
    !is.na(response_pdata$end_time) &
    response_pdata$end_time == "[Completed]"
] <- 0
suppressWarnings(
    response_pdata$end_time <- as.integer(response_pdata$end_time)
)
response_pdata <- response_pdata[order(
    response_pdata$cancer, response_pdata$case_submitter_id,
    response_pdata$drug_name, response_pdata$start_time,
    na.last=TRUE
), ]
response_pdata <- response_pdata[
    !duplicated(response_pdata[c("case_submitter_id", "drug_name")]),
]
response_pdata <- merge(response_pdata, gdc_case_meta, by="case_submitter_id")
response_pdata$cancer <- as.factor(response_pdata$cancer)
response_pdata$response <- as.factor(response_pdata$response)
response_pdata$drug_name <- as.factor(response_pdata$drug_name)

# survival metadata
cat("Processing survival phenotypic data\n")
survival_pdata <- read.delim(
    paste(args$data_dir, "tcga_survival.tsv", sep="/"), row.names=NULL,
    stringsAsFactors=FALSE
)
colnames(survival_pdata)[2] <- "case_submitter_id"
survival_pdata$cancer <- paste("TCGA", survival_pdata$cancer, sep="-")
survival_pdata <- survival_pdata[order(
    survival_pdata$cancer, survival_pdata$case_submitter_id
), ]
survival_pdata <- merge(survival_pdata, gdc_case_meta, by="case_submitter_id")

cat("Writing response_pdata.rds\n")
saveRDS(response_pdata, paste(args$data_dir, "response_pdata.rds", sep="/"))
cat("Writing survival_pdata.rds\n")
saveRDS(survival_pdata, paste(args$data_dir, "survival_pdata.rds", sep="/"))
