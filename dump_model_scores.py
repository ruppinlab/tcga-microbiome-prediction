import os
import re
from argparse import ArgumentParser

import numpy as np
import pandas as pd
import rpy2.rinterface_lib.embedded as r_embedded

r_embedded.set_initoptions(("rpy2", "--quiet", "--no-save", "--max-ppsize=500000"))

import rpy2.robjects as ro
from joblib import dump, load
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

parser = ArgumentParser()
parser.add_argument(
    "--results-dir", type=str, default="results/models", help="results dir"
)
args = parser.parse_args()

metric = {"surv": "score", "resp": "roc_auc"}
surv_model_codes = ["cnet", "cox_clinical"]

r_base = importr("base")

all_scores_dfs = {}
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
            split_results_file = "{}/{}".format(dirpath, filename)
            print("Loading", split_results_file)
            split_results = load(split_results_file)
            scores = []
            for split_result in split_results:
                if split_result is None:
                    scores.append(np.nan)
                else:
                    scores.append(split_result["scores"]["te"][metric[analysis]])
            scores_df = pd.DataFrame({model_name: scores})
            if model_code not in all_scores_dfs:
                all_scores_dfs[model_code] = scores_df
            else:
                all_scores_dfs[model_code] = pd.concat(
                    [all_scores_dfs[model_code], scores_df], axis=1
                )

for model_code, all_scores_df in all_scores_dfs.items():
    analysis = "surv" if model_code in surv_model_codes else "resp"
    out_dir = "{}/{}".format(args.results_dir, analysis)
    os.makedirs(out_dir, mode=0o755, exist_ok=True)
    all_scores_df.to_csv("{}/{}_model_scores.tsv".format(out_dir, model_code), sep="\t")
    dump(all_scores_df, "{}/{}_model_scores.pkl".format(out_dir, model_code))
    with (ro.default_converter + pandas2ri.converter).context():
        r_base.saveRDS(
            all_scores_df, "{}/{}_model_scores.rds".format(out_dir, model_code)
        )
