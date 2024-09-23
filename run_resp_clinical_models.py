import os
import warnings
from argparse import ArgumentParser
from glob import glob
from itertools import product

import numpy as np
import pandas as pd
import rpy2.rinterface_lib.embedded as r_embedded

r_embedded.set_initoptions(("rpy2", "--quiet", "--no-save", "--max-ppsize=500000"))

import rpy2.robjects as ro
from joblib import delayed, dump, Parallel
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr

from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    auc,
    average_precision_score,
    balanced_accuracy_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder, StandardScaler
from sklearn.svm import SVC
from tabulate import tabulate

from sklearn_extensions.model_selection import RepeatedStratifiedGroupKFold


def calculate_test_scores(
    pipe, X_test, y_test, pipe_predict_params, test_sample_weights=None
):
    scores = {}
    if hasattr(pipe, "decision_function"):
        y_score = pipe.decision_function(X_test, **pipe_predict_params)
    else:
        y_score = pipe.predict_proba(X_test, **pipe_predict_params)[:, 1]
    scores["y_score"] = y_score
    for metric in metrics:
        if metric == "roc_auc":
            scores[metric] = roc_auc_score(
                y_test, y_score, sample_weight=test_sample_weights
            )
            scores["fpr"], scores["tpr"], _ = roc_curve(
                y_test, y_score, pos_label=1, sample_weight=test_sample_weights
            )
        elif metric == "balanced_accuracy":
            y_pred = pipe.predict(X_test, **pipe_predict_params)
            scores["y_pred"] = y_pred
            scores[metric] = balanced_accuracy_score(
                y_test, y_pred, sample_weight=test_sample_weights
            )
        elif metric == "average_precision":
            scores[metric] = average_precision_score(
                y_test, y_score, sample_weight=test_sample_weights
            )
            scores["pre"], scores["rec"], _ = precision_recall_curve(
                y_test, y_score, pos_label=1, sample_weight=test_sample_weights
            )
            scores["pr_auc"] = auc(scores["rec"], scores["pre"])
    return scores


def fit_models(
    pipe, X, y, groups, sample_weights, test_splits, test_repeats, dataset_name
):
    if args.n_jobs == 1 and args.verbose > 1:
        model_code = "svm" if isinstance(pipe[-1], SVC) else "lgr"
        model_name = "_".join([dataset_name, model_code, "clinical"])
        print(model_name)

    if groups is None:
        cv = RepeatedStratifiedKFold(
            n_splits=test_splits, n_repeats=test_repeats, random_state=random_seed
        )
    else:
        cv = RepeatedStratifiedGroupKFold(
            n_splits=test_splits, n_repeats=test_repeats, random_state=random_seed
        )

    split_results = []
    for train_idxs, test_idxs in cv.split(X, y, groups):
        train_sample_weights = None
        test_sample_weights = None
        if sample_weights is not None:
            train_sample_weights = sample_weights[train_idxs]
            test_sample_weights = sample_weights[test_idxs]
        pipe.fit(
            X.iloc[train_idxs], y[train_idxs], clf1__sample_weight=train_sample_weights
        )
        split_scores = {
            "te": calculate_test_scores(
                pipe,
                X.iloc[test_idxs],
                y[test_idxs],
                pipe_predict_params={},
                test_sample_weights=test_sample_weights,
            )
        }
        split_results.append({"scores": split_scores})

    return split_results


parser = ArgumentParser()
parser.add_argument("--data-dir", type=str, default="data", help="data dir")
parser.add_argument(
    "--results-dir", type=str, default="results/models", help="results dir"
)
parser.add_argument("--test-splits", type=int, help="num test splits")
parser.add_argument("--test-repeats", type=int, help="num test repeats")
parser.add_argument("--n-jobs", type=int, default=-1, help="num parallel jobs")
parser.add_argument(
    "--parallel-backend", type=str, default="loky", help="joblib parallel backend"
)
parser.add_argument(
    "--show-warnings", default=False, action="store_true", help="show fit warnings"
)
parser.add_argument("--verbose", type=int, default=1, help="verbosity")
args = parser.parse_args()

random_seed = 777

if not args.show_warnings:
    if args.parallel_backend == "multiprocessing":
        warnings.filterwarnings(
            "ignore",
            category=ConvergenceWarning,
            message=(
                "^The max_iter was reached which means the coef_ did not " "converge"
            ),
            module="sklearn.linear_model._sag",
        )
    else:
        python_warnings = (
            [os.environ["PYTHONWARNINGS"]] if "PYTHONWARNINGS" in os.environ else []
        )
        python_warnings.append(
            ":".join(
                [
                    "ignore",
                    (
                        "The max_iter was reached which means the coef_ did not "
                        "converge"
                    ),
                    "UserWarning",
                    "sklearn.linear_model._sag",
                ]
            )
        )
        os.environ["PYTHONWARNINGS"] = ",".join(python_warnings)

out_dir = "{}/resp".format(args.results_dir)
os.makedirs(out_dir, mode=0o755, exist_ok=True)

r_base = importr("base")
r_biobase = importr("Biobase")

metrics = ["roc_auc", "average_precision", "balanced_accuracy"]
ordinal_encoder_categories = {
    "tumor_stage": ["i", "ii", "iii", "iv"],
}

pipes = [
    Pipeline(
        [
            ("trf0", StandardScaler()),
            (
                "clf1",
                SVC(
                    kernel="linear",
                    class_weight="balanced",
                    random_state=random_seed,
                ),
            ),
        ]
    ),
    Pipeline(
        [
            ("trf0", StandardScaler()),
            (
                "clf1",
                LogisticRegression(
                    penalty="l2",
                    solver="lbfgs",
                    max_iter=5000,
                    class_weight="balanced",
                    random_state=random_seed,
                ),
            ),
        ]
    ),
]

datasets = []
eset_files = sorted(glob("{}/tcga_*_resp_*_eset.rds".format(args.data_dir)))
num_esets = len(eset_files)
for eset_idx, eset_file in enumerate(eset_files):
    file_basename = os.path.splitext(os.path.split(eset_file)[1])[0]
    _, cancer, _, target, _, *rest = file_basename.split("_")

    cancer_target = "_".join([cancer, target])
    if args.test_splits is None:
        test_splits = (
            3
            if cancer_target
            in [
                "blca_doxorubicin",
                "brca_anastrazole",
                "brca_tamoxifen",
                "coad_bevacizumab",
                "coad_irinotecan",
                "esca_capecitabine",
                "lusc_cisplatin",
                "lusc_docetaxel",
                "lusc_gemcitabine",
                "sarc_ifosfamide",
            ]
            else 4
        )
    else:
        test_splits = args.test_splits
    if args.test_repeats is None:
        test_repeats = (
            33
            if cancer_target
            in [
                "blca_doxorubicin",
                "brca_anastrazole",
                "brca_tamoxifen",
                "coad_bevacizumab",
                "coad_irinotecan",
                "esca_capecitabine",
                "lusc_cisplatin",
                "lusc_docetaxel",
                "lusc_gemcitabine",
                "sarc_ifosfamide",
            ]
            else 25
        )
    else:
        test_repeats = args.test_repeats

    dataset_name = "_".join(file_basename.split("_")[:-1])

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
    y = np.array(sample_meta["Class"], dtype=int)

    if "Group" in sample_meta.columns:
        groups = np.array(sample_meta["Group"], dtype=int)
        _, group_indices, group_counts = np.unique(
            groups, return_inverse=True, return_counts=True
        )
        sample_weights = (np.max(group_counts) / group_counts)[group_indices]
    else:
        groups = None
        sample_weights = None

    X["age_at_diagnosis"] = sample_meta[["age_at_diagnosis"]]
    if sample_meta["gender"].nunique() > 1:
        ohe = OneHotEncoder(drop="first", sparse=False)
        ohe.fit(sample_meta[["gender"]])
        feature_name = "gender_{}".format(ohe.categories_[0][1])
        X[feature_name] = ohe.transform(sample_meta[["gender"]])
    if sample_meta["tumor_stage"].nunique() > 1:
        ode = OrdinalEncoder(categories=[ordinal_encoder_categories["tumor_stage"]])
        ode.fit(sample_meta[["tumor_stage"]])
        X["tumor_stage"] = ode.transform(sample_meta[["tumor_stage"]])

    datasets.append(
        (X, y, groups, sample_weights, test_splits, test_repeats, dataset_name)
    )

if args.verbose < 2:
    print(flush=True)

print("Running drug response clinical models", flush=True)
all_results = Parallel(
    n_jobs=args.n_jobs, backend=args.parallel_backend, verbose=args.verbose
)(
    delayed(fit_models)(
        pipe, X, y, groups, sample_weights, test_splits, test_repeats, dataset_name
    )
    for pipe, (
        X,
        y,
        groups,
        sample_weights,
        test_splits,
        test_repeats,
        dataset_name,
    ) in (product(pipes, datasets))
)

if args.verbose < 1:
    print(flush=True)

mean_scores = []
for (pipe, eset_file), split_results in zip(product(pipes, eset_files), all_results):
    file_basename = os.path.splitext(os.path.split(eset_file)[1])[0]
    _, cancer, analysis, target, data_type, *rest = file_basename.split("_")

    roc_scores, pr_scores = [], []
    for split_result in split_results:
        roc_scores.append(split_result["scores"]["te"]["roc_auc"])
        pr_scores.append(split_result["scores"]["te"]["pr_auc"])

    dataset_name = "_".join(file_basename.split("_")[:-1])
    model_code = "svm" if isinstance(pipe[-1], SVC) else "lgr"
    model_name = "_".join([dataset_name, model_code, "clinical"])

    mean_score = np.nanmean(roc_scores)
    mean_scores.append([analysis, cancer, target, data_type, model_code, mean_score])

    results_dir = "{}/{}".format(out_dir, model_name)
    os.makedirs(results_dir, mode=0o755, exist_ok=True)
    dump(split_results, "{}/{}_split_results.pkl".format(results_dir, model_name))

mean_scores_df = pd.DataFrame(
    mean_scores,
    columns=["Analysis", "Cancer", "Target", "Data Type", "Model Code", "Mean Score"],
)
mean_scores_df.to_csv(
    "{}/resp_clinical_model_mean_scores.tsv".format(args.results_dir),
    index=False,
    sep="\t",
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
