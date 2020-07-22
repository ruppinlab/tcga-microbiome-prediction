# Authors: Leandro Cruz Hermida <hermidal@cs.umd.edu>

suppressPackageStartupMessages(library("edgeR"))


edger_filterbyexpr_mask <- function(
    X, y, sample_meta=NULL, model_batch=FALSE, is_classif=TRUE
) {
    dge <- DGEList(counts=t(X))
    if (
        model_batch && !is.null(sample_meta) &&
        length(unique(sample_meta$Batch)) > 1
    ) {
        sample_meta$Batch <- factor(sample_meta$Batch)
        if (is_classif) {
            sample_meta$Class <- factor(sample_meta$Class)
            design <- model.matrix(~Batch + Class, data=sample_meta)
        } else {
            design <- model.matrix(~Batch, data=sample_meta)
        }
    } else if (is_classif) {
        design <- model.matrix(~factor(y))
    } else {
        design <- NULL
    }
    return(filterByExpr(dge, design))
}
