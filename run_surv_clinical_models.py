import os
import warnings
from argparse import ArgumentParser
from glob import glob

import numpy as np
import pandas as pd
import rpy2.rinterface_lib.embedded as r_embedded

r_embedded.set_initoptions(("rpy2", "--quiet", "--no-save", "--max-ppsize=500000"))

import rpy2.robjects as ro
from joblib import delayed, dump, Parallel
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder, StandardScaler
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.util import Surv
from tabulate import tabulate

from sksurv_extensions.model_selection import (
    SurvivalStratifiedShuffleSplit,
    SurvivalStratifiedSampleFromGroupShuffleSplit,
)


def fit_models(X, y, groups, group_weights, test_splits, test_size):
    pipe = Pipeline(
        [
            ("trf0", StandardScaler()),
            (
                "srv1",
                CoxPHSurvivalAnalysis(
                    alpha=1e-09, n_iter=10000, ties="efron", tol=1e-09
                ),
            ),
        ]
    )

    if groups is None:
        cv = SurvivalStratifiedShuffleSplit(
            n_splits=test_splits, test_size=test_size, random_state=random_seed
        )
        cv_split_params = {}
    else:
        cv = SurvivalStratifiedSampleFromGroupShuffleSplit(
            n_splits=test_splits, test_size=test_size, random_state=random_seed
        )
        cv_split_params = {"weights": group_weights}

    split_models, split_results = [], []
    for train_idxs, test_idxs in cv.split(X, y, groups, **cv_split_params):
        try:
            pipe.fit(X.iloc[train_idxs], y[train_idxs])
            score = pipe.score(X.iloc[test_idxs], y[test_idxs])
            y_pred = pipe.predict(X.iloc[test_idxs])
            surv_funcs = pipe.predict_survival_function(X.iloc[test_idxs])
        except Exception:
            split_models.append(None)
            split_results.append(None)
        else:
            split_models.append(pipe)
            split_scores = {"te": {"score": score, "y_pred": y_pred}}
            split_results.append({"scores": split_scores, "surv_funcs": surv_funcs})

    return split_models, split_results


parser = ArgumentParser()
parser.add_argument("--data-dir", type=str, default="data", help="data dir")
parser.add_argument(
    "--results-dir", type=str, default="results/models", help="results dir"
)
parser.add_argument("--test-splits", type=int, help="num test splits")
parser.add_argument("--test-size", type=float, help="test split size")
parser.add_argument("--n-jobs", type=int, default=-1, help="num parallel jobs")
parser.add_argument(
    "--parallel-backend", type=str, default="loky", help="joblib parallel backend"
)
parser.add_argument("--verbose", type=int, default=1, help="verbosity")
args = parser.parse_args()

test_splits = 100 if args.test_splits is None else args.test_splits
test_size = 0.25 if args.test_size is None else args.test_size
random_seed = 777

if args.parallel_backend == "multiprocessing":
    warnings.filterwarnings(
        "ignore",
        category=RuntimeWarning,
        message="^divide by zero encountered in true_divide",
        module="sksurv.linear_model.coxph",
    )
    warnings.filterwarnings(
        "ignore",
        category=RuntimeWarning,
        message="^invalid value encountered in divide",
    )
    warnings.filterwarnings(
        "ignore",
        category=RuntimeWarning,
        message="^overflow encountered in exp",
        module="sksurv.linear_model.coxph",
    )
    warnings.filterwarnings(
        "ignore",
        category=RuntimeWarning,
        message="^overflow encountered in power",
        module="sksurv.linear_model.coxph",
    )
else:
    python_warnings = (
        [os.environ["PYTHONWARNINGS"]] if "PYTHONWARNINGS" in os.environ else []
    )
    python_warnings.append(
        ":".join(
            [
                "ignore",
                "divide by zero encountered in true_divide",
                "RuntimeWarning",
                "sksurv.linear_model.coxph",
            ]
        )
    )
    python_warnings.append(
        ":".join(
            [
                "ignore",
                "invalid value encountered in divide",
                "RuntimeWarning",
            ]
        )
    )
    python_warnings.append(
        ":".join(
            [
                "ignore",
                "overflow encountered in exp",
                "RuntimeWarning",
                "sksurv.linear_model.coxph",
            ]
        )
    )
    python_warnings.append(
        ":".join(
            [
                "ignore",
                "overflow encountered in power",
                "RuntimeWarning",
                "sksurv.linear_model.coxph",
            ]
        )
    )
    os.environ["PYTHONWARNINGS"] = ",".join(python_warnings)

out_dir = "{}/surv".format(args.results_dir)
os.makedirs(out_dir, mode=0o755, exist_ok=True)

ordinal_encoder_categories = {
    "tumor_stage": ["i", "ii", "iii", "iv"],
}

r_base = importr("base")
r_biobase = importr("Biobase")

datasets = []
eset_files = sorted(glob("{}/tcga_*_surv_*_eset.rds".format(args.data_dir)))
num_esets = len(eset_files)
for eset_idx, eset_file in enumerate(eset_files):
    file_basename = os.path.splitext(os.path.split(eset_file)[1])[0]

    if args.verbose < 2:
        print(
            "Loading {:d}/{:d} esets".format(eset_idx + 1, num_esets),
            end="\r",
            flush=True,
        )
    else:
        print("Loading {}".format(file_basename))

    eset = r_base.readRDS(eset_file)
    with (ro.default_converter + numpy2ri.converter + pandas2ri.converter).context():
        sample_meta = r_biobase.pData(eset)
        X = pd.DataFrame(index=sample_meta.index)
        y = Surv.from_dataframe("Status", "Survival_in_days", sample_meta)

    if "Group" in sample_meta.columns:
        groups = np.array(sample_meta["Group"], dtype=int)
        if (
            "GroupWeight" in sample_meta.columns
            and sample_meta["GroupWeight"].unique().size > 1
        ):
            group_weights = np.array(sample_meta["GroupWeight"], dtype=float)
        else:
            group_weights = None
    else:
        groups = None
        group_weights = None

    X["age_at_diagnosis"] = sample_meta[["age_at_diagnosis"]]
    if sample_meta["gender"].unique().size > 1:
        ohe = OneHotEncoder(drop="first", sparse=False)
        ohe.fit(sample_meta[["gender"]])
        feature_name = "gender_{}".format(ohe.categories_[0][1])
        X[feature_name] = ohe.transform(sample_meta[["gender"]])
    if sample_meta["tumor_stage"].unique().size > 1:
        ode = OrdinalEncoder(categories=[ordinal_encoder_categories["tumor_stage"]])
        ode.fit(sample_meta[["tumor_stage"]])
        X["tumor_stage"] = ode.transform(sample_meta[["tumor_stage"]])

    datasets.append((X, y, groups, group_weights))

if args.verbose < 2:
    print(flush=True)

print("Running survival clinical models", flush=True)
all_models, all_results = zip(
    *Parallel(n_jobs=args.n_jobs, backend=args.parallel_backend, verbose=args.verbose)(
        delayed(fit_models)(X, y, groups, group_weights, test_splits, test_size)
        for X, y, groups, group_weights in datasets
    )
)

if args.verbose < 1:
    print(flush=True)

mean_scores = []
for eset_file, split_models, split_results in zip(eset_files, all_models, all_results):
    file_basename = os.path.splitext(os.path.split(eset_file)[1])[0]
    _, cancer, analysis, target, data_type, *rest = file_basename.split("_")

    scores = []
    for split_result in split_results:
        if split_result is not None:
            scores.append(split_result["scores"]["te"]["score"])
        else:
            scores.append(np.nan)

    dataset_name = "_".join(file_basename.split("_")[:-1])
    model_code = "cox"
    model_name = "_".join([dataset_name, model_code, "clinical"])

    mean_score = np.nanmean(scores)
    mean_scores.append([analysis, cancer, target, data_type, model_code, mean_score])

    results_dir = "{}/{}".format(out_dir, model_name)
    os.makedirs(results_dir, mode=0o755, exist_ok=True)
    dump(split_models, "{}/{}_split_models.pkl".format(results_dir, model_name))
    dump(split_results, "{}/{}_split_results.pkl".format(results_dir, model_name))

mean_scores_df = pd.DataFrame(
    mean_scores,
    columns=["Analysis", "Cancer", "Target", "Data Type", "Model Code", "Mean Score"],
)
mean_scores_df.to_csv(
    "{}/surv_clinical_model_mean_scores.tsv".format(out_dir), index=False, sep="\t"
)
if args.verbose > 0:
    print(
        tabulate(
            mean_scores_df.sort_values(
                ["Analysis", "Cancer", "Target", "Data Type", "Model Code"]
            ),
            floatfmt=".4f",
            showindex=False,
            headers="keys",
        )
    )
