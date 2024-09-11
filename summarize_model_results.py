import os
import re
from argparse import ArgumentParser
from glob import glob

import numpy as np
import pandas as pd
from joblib import load

parser = ArgumentParser()
parser.add_argument(
    "--results-dir", type=str, default="results/models", help="results dir"
)
args = parser.parse_args()

surv_mean_scores = "{}/surv/surv_clinical_model_mean_scores.tsv".format(
    args.results_dir
)
resp_mean_scores = "{}/resp/resp_clinical_model_mean_scores.tsv".format(
    args.results_dir
)

metric = {"surv": "score", "resp": "roc_auc"}
penalty_factor_meta_col = "Penalty Factor"
lgr_model_codes = ["edge", "elgr", "voom", "zinb"]

results = []
split_results_regex = re.compile("^(.+?)_split_results\\.pkl$")
for dirpath, dirnames, filenames in sorted(os.walk(args.results_dir)):
    for filename in filenames:
        if m := re.search(split_results_regex, filename):
            model_name = m.group(1)
            _, cancer, analysis, target, data_type, *rest = model_name.split("_")
            if data_type == "htseq":
                model_code = "_".join(rest[1:])
            else:
                model_code = "_".join(rest)

            if model_code.endswith("clinical"):
                continue

            job_id = ""
            if s := glob("{}/slurm-*.out".format(dirpath)):
                slurm_basename = os.path.splitext(os.path.split(s[0])[1])[0]
                job_id = re.findall("\\d+", slurm_basename)[0]

            split_results_file = "{}/{}".format(dirpath, filename)
            print("Loading", split_results_file)
            split_results = load(split_results_file)

            scores = []
            num_features = []
            for split_result in split_results:
                if split_result is None:
                    continue
                scores.append(split_result["scores"]["te"][metric[analysis]])
                split_feature_meta = split_result["feature_meta"]
                if penalty_factor_meta_col in split_feature_meta.columns:
                    num_features.append(
                        split_feature_meta.loc[
                            split_feature_meta[penalty_factor_meta_col] != 0
                        ].shape[0]
                    )
                else:
                    num_features.append(split_feature_meta.shape[0])

            results.append(
                [
                    analysis,
                    cancer,
                    target,
                    data_type,
                    model_code,
                    np.mean(scores),
                    np.mean(num_features),
                    job_id,
                ]
            )

results_summary = pd.DataFrame(
    results,
    columns=[
        "Analysis",
        "Cancer",
        "Target",
        "Data Type",
        "Model Code",
        "Mean Score",
        "Mean Num Features",
        "Job ID",
    ],
)
if os.path.isfile(surv_mean_scores):
    surv_mean_scores = pd.read_csv(surv_mean_scores, sep="\t")
    surv_mean_scores["Model Code"] = surv_mean_scores["Model Code"].replace(
        "cox", "cnet"
    )
    results_summary = pd.merge(
        results_summary,
        surv_mean_scores,
        how="left",
        on=["Analysis", "Cancer", "Target", "Data Type", "Model Code"],
    )
if os.path.isfile(resp_mean_scores):
    resp_mean_scores = pd.read_csv(resp_mean_scores, sep="\t")
    resp_mean_scores["Model Code"] = resp_mean_scores["Model Code"].replace(
        "svm", "srfe"
    )
    lgr_mean_scores = resp_mean_scores.loc[
        resp_mean_scores["Model Code"] == "lgr"
    ].copy()
    for model_code in lgr_model_codes:
        lgr_mean_scores["Model Code"] = model_code
        resp_mean_scores = pd.concat([resp_mean_scores, lgr_mean_scores], axis=0)
    results_summary = pd.merge(
        results_summary,
        resp_mean_scores,
        how="left",
        on=["Analysis", "Cancer", "Target", "Data Type", "Model Code"],
    )
results_summary["Clinical Mean Score"] = results_summary["Mean Score_y"].combine_first(
    results_summary["Mean Score"]
)
results_summary.drop(columns=["Mean Score", "Mean Score_y"], inplace=True)
results_summary.rename(columns={"Mean Score_x": "Mean Score"}, inplace=True)
results_summary.sort_values(
    by=["Analysis", "Cancer", "Target", "Data Type", "Model Code", "Mean Score"],
    inplace=True,
)
results_summary["Mean Num Features"] = results_summary["Mean Num Features"].round()
results_summary.to_csv(
    "{}/model_results_summary.tsv".format(args.results_dir),
    sep="\t",
    index=False,
    float_format="%.4f",
)
