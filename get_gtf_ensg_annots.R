suppressPackageStartupMessages({
    library(AnnotationHub)
    library(argparser)
    library(ensembldb)
    library(httr)
    library(rtracklayer)
    library(stringr)
    library(tools)
})

argp <- arg_parser("Get GDC GENCODE gtf and latest Ensembl annotations")
argp <- add_argument(
    argp, "--data-dir",
    type = "character", default = "data", help = "Data directory"
)
argp <- add_argument(
    argp, "--ensdb-version",
    default = 102, help = "Ensembl release version"
)
args <- parse_args(argp)

gdc_data_base_url <- "https://api.gdc.cancer.gov/data/"
gdc_gtf_url <- paste0(gdc_data_base_url, "be002a2c-3b27-43f3-9e0f-fd47db92a6b5")
gdc_gtf_filename <- "gencode.v36.annotation.gtf.gz"
gdc_gtf_file <- paste(args$data_dir, gdc_gtf_filename, sep = "/")
if (file.exists(gdc_gtf_file)) {
    cat("Using existing", gdc_gtf_file, "\n")
} else {
    cat("Downloading", gdc_gtf_file, "\n")
    download.file(gdc_gtf_url, gdc_gtf_file)
}
if (md5sum(gdc_gtf_file) != "c03931958d4572148650d62eb6dec41a") {
    stop(gdc_gtf_filename, "md5 doesn't match")
}

cat("Loading", basename(gdc_gtf_file), "\n")
gencode_gtf <- import(gdc_gtf_file)
gene_ids <- sort(unique(gencode_gtf$gene_id))
gene_ids_no_version <- unlist(
    strsplit(gene_ids, ".", fixed = TRUE)
)[2 * (seq_len(length(gene_ids))) - 1]

ah_query_pattern <- c("Homo sapiens", "EnsDb", args$ensdb_version)
cat("Getting", ah_query_pattern, "\n")
suppressWarnings(ah <- AnnotationHub())
ahdb <- query(ah, pattern = ah_query_pattern)
edb <- ahdb[[1]]

cat("Processing annotations", "\n")
edb_gene_id_annots <- select(
    edb,
    keys = gene_ids_no_version, keytype = "GENEID",
    columns = c("GENEID", "SYMBOL"), order.by = "GENEID"
)
edb_gene_id_annots <- data.frame(
    Symbol = edb_gene_id_annots$SYMBOL, row.names = edb_gene_id_annots$GENEID
)
gene_id_symbols <- sapply(
    gene_ids_no_version,
    function(gene_id) {
        if (gene_id %in% row.names(edb_gene_id_annots)) {
            as.character(edb_gene_id_annots[gene_id, "Symbol"])
        } else {
            ""
        }
    }
)

out_filename <- paste(
    "gencode", str_match(gdc_gtf_filename, "\\.(v\\d+)\\.")[1, 2],
    "ensg", paste0("v", args$ensdb_version), "annots.tsv",
    sep = "_"
)
out_file <- paste(args$data_dir, out_filename, sep = "/")
cat("Writing", out_file, "\n")
write.table(
    data.frame(ID_REF = gene_ids, Symbol = gene_id_symbols),
    file = out_file, quote = FALSE, row.names = FALSE, sep = "\t"
)
