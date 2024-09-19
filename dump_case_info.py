import os
from argparse import ArgumentParser
from glob import glob

import numpy as np
import pandas as pd
import rpy2.rinterface_lib.embedded as r_embedded

r_embedded.set_initoptions(("rpy2", "--quiet", "--no-save", "--max-ppsize=500000"))

import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr
from tabulate import tabulate

parser = ArgumentParser()
parser.add_argument("--data-dir", type=str, default="data", help="data dir")
parser.add_argument(
    "--results-dir", type=str, default="results/models", help="results dir"
)
parser.add_argument(
    "--resp-sort-by",
    type=str,
    nargs="+",
    choices=[
        "Cancer",
        "Analysis",
        "Target",
        "Data Type",
        "Num Cases",
        "NR (-) Cases",
        "R (+) Cases",
    ],
    default=["Cancer", "Analysis", "Target", "Data Type"],
    help="Drug response columns to sort by",
)
parser.add_argument(
    "--surv-sort-by",
    type=str,
    nargs="+",
    choices=["Cancer", "Analysis", "Target", "Data Type", "Num Cases"],
    default=["Cancer", "Analysis", "Target", "Data Type"],
    help="Survival columns to sort by",
)
args = parser.parse_args()

os.makedirs(args.results_dir, mode=0o755, exist_ok=True)

r_base = importr("base")
r_biobase = importr("Biobase")

results = []
eset_files = sorted(glob("{}/tcga_*_resp_*_eset.rds".format(args.data_dir)))
num_esets = len(eset_files)
for eset_idx, eset_file in enumerate(eset_files):
    file_basename = os.path.splitext(os.path.split(eset_file)[1])[0]
    _, cancer, analysis, target, data_type, *rest = file_basename.split("_")
    print(
        "Loading {:d}/{:d} esets".format(eset_idx + 1, num_esets), end="\r", flush=True
    )
    eset = r_base.readRDS(eset_file)
    with (ro.default_converter + numpy2ri.converter + pandas2ri.converter).context():
        sample_meta = r_biobase.pData(eset)
    num_cases = sample_meta["case_submitter_id"].nunique()
    neg_cases, pos_cases = sample_meta.groupby("Class", observed=False)[
        "case_submitter_id"
    ].nunique()
    cancer = cancer.upper()
    target = target.title()
    data_type = (
        "Microbiome"
        if data_type == "kraken"
        else "Expression" if data_type == "htseq" else "Combo"
    )
    results.append(
        [cancer, analysis, target, data_type, num_cases, neg_cases, pos_cases]
    )

results_df = pd.DataFrame(
    results,
    columns=[
        "Cancer",
        "Analysis",
        "Target",
        "Data Type",
        "Num Cases",
        "NR (-) Cases",
        "R (+) Cases",
    ],
)
results_df.to_csv(
    "{}/resp_case_info.tsv".format(args.results_dir), sep="\t", index=False
)
print(
    tabulate(
        results_df.sort_values(by=args.resp_sort_by), headers="keys", showindex=False
    )
)

print()

results = []
eset_files = sorted(glob("{}/tcga_*_surv_*_eset.rds".format(args.data_dir)))
num_esets = len(eset_files)
for eset_idx, eset_file in enumerate(eset_files):
    file_basename = os.path.splitext(os.path.split(eset_file)[1])[0]
    _, cancer, analysis, target, data_type, *rest = file_basename.split("_")
    print(
        "Loading {:d}/{:d} esets".format(eset_idx + 1, num_esets), end="\r", flush=True
    )
    eset = r_base.readRDS(eset_file)
    with (ro.default_converter + numpy2ri.converter + pandas2ri.converter).context():
        sample_meta = r_biobase.pData(eset)
    num_cases = sample_meta["case_submitter_id"].nunique()
    cancer = cancer.upper()
    target = target.upper()
    data_type = (
        "Microbiome"
        if data_type == "kraken"
        else "Expression" if data_type == "htseq" else "Combo"
    )
    results.append([cancer, analysis, target, data_type, num_cases])

results_df = pd.DataFrame(
    results, columns=["Cancer", "Analysis", "Target", "Data Type", "Num Cases"]
)
results_df.to_csv(
    "{}/surv_case_info.tsv".format(args.results_dir), sep="\t", index=False
)
print(
    tabulate(
        results_df.sort_values(by=args.surv_sort_by), headers="keys", showindex=False
    )
)
