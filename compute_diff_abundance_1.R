suppressPackageStartupMessages({
    library(Biobase)
    library(edgeR)
    library(EnhancedVolcano)
    library(stringr)
    library(Wrench)
    library(zinbwave)
})

fc <- 1.0
lfc <- log2(fc)
padj <- 0.05
padj_meth <- "BH"
fig_dim <- 8
fig_dpi <- 300

for (eset_file in Sys.glob("data/tcga_*_resp_*_kraken_eset.rds")) {
    cat("Loading", eset_file, "\n")
    eset <- readRDS(eset_file)

    counts <- exprs(eset)
    pdata <- pData(eset)
    fdata <- fData(eset)

    dge <- DGEList(counts = counts, samples = pdata, genes = fdata)
    dge <- sumTechReps(dge, ID = dge$samples$case_submitter_id)

    # W <- wrench(dge$counts, condition = dge$samples$Class)
    # dge$samples$norm.factors <- W$ccf

    dge$samples$Class <- as.factor(dge$samples$Class)
    design_formula <- ~Class
    design <- model.matrix(design_formula, data = dge$samples)

    # dge <- calcNormFactors(dge, method = "TMMwsp")
    # dge <- estimateDisp(dge, design, robust = TRUE)
    # fit <- glmQLFit(dge, design, robust = TRUE)
    # glt <- glmTreat(fit, coef = ncol(design), lfc = lfc)
    # results <- as.data.frame(
    #     topTags(glt, n = Inf, adjust.method = padj_meth, sort.by = "PValue")
    # )

    zinb <- zinbwave(
        SummarizedExperiment(
            assays = list(counts = dge$counts), colData = dge$samples
        ),
        X = design_formula, K = 0, epsilon = 1e12, zeroinflation = TRUE,
        observationalWeights = TRUE
    )
    dge <- calcNormFactors(dge, method = "TMM")
    dge$weights <- assay(zinb, "weights")
    dge <- estimateDisp(dge, design, robust = TRUE)
    fit <- glmFit(dge, design)
    lrt <- glmWeightedF(fit, coef = ncol(design), independentFiltering = FALSE)
    results <- as.data.frame(
        topTags(lrt, n = Inf, adjust.method = padj_meth, sort.by = "PValue")
    )
    # results$PValue[is.na(results$PValue)] <- 1
    # results$padjFilter[is.na(results$padjFilter)] <- 1
    # results$FDR[is.na(results$FDR)] <- 1

    filename_parts <- str_split(basename(eset_file), "_")[[1]]
    cancer <- filename_parts[2]
    analysis <- filename_parts[3]
    target <- filename_parts[4]

    p <- EnhancedVolcano(
        results,
        lab = row.names(results),
        x = "logFC",
        y = "PValue",
        xlim = c(floor(min(results$logFC)), ceiling(max(results$logFC))),
        ylim = c(0, ceiling(max(-log10(results$PValue)))),
        pCutoff = (
            tail(results$PValue[results$FDR < padj], n = 1)
            + head(results$PValue[results$FDR > padj], n = 1)
        ) / 2,
        FCcutoff = lfc,
        pointSize = 3.0,
        labSize = 3.5,
        labFace = "bold",
        drawConnectors = TRUE,
        widthConnectors = 1,
        colConnectors = "grey20",
        maxoverlapsConnectors = 20,
        arrowheads = FALSE,
        boxedLabels = TRUE,
        title = paste(str_to_upper(cancer), target),
        subtitle = NULL,
        caption = NULL,
        legendPosition = "none"
    )

    dir_path <- "figures/volcano"
    if (!dir.exists(dir_path)) {
        dir.create(
            dir_path,
            showWarnings = FALSE, recursive = TRUE, mode = "0755"
        )
    }
    ggsave(
        file = paste(
            dir_path,
            paste("tcga", cancer, analysis, target, "volcano.png", sep = "_"),
            sep = "/"
        ),
        plot = p, device = "png", width = fig_dim,
        height = fig_dim, units = "in", dpi = fig_dpi
    )
}
