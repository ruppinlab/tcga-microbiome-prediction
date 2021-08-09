# results/analysis

The results/analysis directory contains the results of applying
statistical tests to the ML runs, choosing successful runs and
successful cancer/target pairs, and choosing features.  The files of
most interest are

    - goodness_hits.txt
    - selected_hits.txt
    - simplified_microbial.txt
    - simplified_expression.txt
    - microbial_feature_stats.txt

### goodness_hits.txt
The file goodness_hits.txt is a table of models of both survival and
drug response that meet our criteria to be considered hits.  The
measure of "goodness" is C-index for survival runs and AUROC for drug
response runs, and because the file contains both, it was given the
vague name "goodness_hits.txt".  Columns in the table are

    - cancer - TCGA cancer code
    - analysis - one of resp = Response, surv = Survival
    - versus - OS or PFI for survival, drug name for response
    - features - kraken=microbial, htseq=expression, combo=both
    - how - Machine learning method
    - avg_test - average C-index or AUROC
    - sd_test - standard deviation of C-index or AUROC
    - avg_cov - average C-index or AUROC of the clinical cov model
    - sd_cov  - average C-index or AUROC of the clinical cov model
    - p_value - p-value of a two-sided Wilcoxon signed rank test
    - p_greater - p-value of a Wilcoxon test that the clinical model is worse
    - p_adj - FDR corrected p_value

Though the full ML models are usually better than the clinical
covariate only models, this is false often enough that we did not feel
that using the one-sided Wilcoxon test was justified.  However,
p_greater can correctly determine, conditional on one model being
better, which model is better.

### selected_hits.txt
A simplified version of goodness_hits.txt that also applies the rule
that for drug response 2/3 hits must agree.  So not all models in
goodness_hits.txt appear in selected_hits.txt.

### simplified_microbial.txt
A table of microbial features, listed by cancer/target, for which at
least two ML methods listed in selected_hits.txt call the feature.
Also shows univariate analysis of the microbial features.

## simplified_expression.txt
A table of expression features, listed by cancer/target, for which at
least two ML methods listed in selected_hits.txt call the feature.

### microbial_feature_stats.txt
A relatively free-form file that shows various summaries and
statistics for the features in simplified_microbial.txt.

## Intermediate files

The additional files in the directory are intermediate files.  The
intermediate files that may be of some interest are

    - compared_runs.txt
    - potential_hits.txt
    - expression_features.txt
    - microbial_features.txt

### compared_runs.txt
All runs compared without filtering for quality.

### potential_hits.txt
Runs filtered so that the measure of goodness -- AUROC or C-index --
is at least 0.6 and so that the inference is in the "right" direction,
i.e.  so that the full ML model is better than the covariate-only model.
All filters have been applied except for the filter on FDR.  This file
is of some use to see which models have good AUROC or C-index, but are
excluded because the full models are not consistently better than the
clinical covariate only models.

### expression_features.txt
Expression features before the 2/3 rule is applied.  The raw data for
`simplified_expression.txt`.

### microbial_features.txt
Microbial features before the 2/3 rule is applied.  The raw data for
`simplified_microbial.txt`.
