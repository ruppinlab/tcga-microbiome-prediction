import os
import re
import sys
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from joblib import load
from matplotlib import ticker
from matplotlib.offsetbox import AnchoredText
from scipy.stats import iqr

# suppress linux conda qt5 wayland warning
if sys.platform.startswith("linux"):
    os.environ["XDG_SESSION_TYPE"] = "x11"

parser = ArgumentParser()
parser.add_argument("--results-dir", type=str, default="results", help="results dir")
parser.add_argument("--out-dir", type=str, default="figures/perm_hist", help="out dir")
parser.add_argument(
    "--model-code",
    type=str,
    nargs="+",
    choices=["edge", "elgr", "srfe", "voom", "zinb"],
    default=["edge", "elgr", "srfe", "voom", "zinb"],
    help="response model code filter",
)
parser.add_argument(
    "--file-format",
    type=str,
    nargs="+",
    choices=["png", "pdf", "svg", "tif"],
    default=["pdf"],
    help="save file format",
)
args = parser.parse_args()

model_results_dir = "{}/models".format(args.results_dir)

os.makedirs(args.out_dir, mode=0o755, exist_ok=True)

title_fontsize = 22
axis_fontsize = 20
legend_fontsize = 20
fig_dim = 4
fig_dpi = 300

plt.rcParams["figure.max_open_warning"] = 0
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "Helvetica",
    "Nimbus Sans",
    "Arial",
    "DejaVu Sans",
    "sans-serif",
]

model_codes_regex = "|".join(args.model_code)
perm_results_regex = re.compile(
    "^(.+?_(?:{}))_perm_results\\.pkl$".format(model_codes_regex)
)
for dirpath, dirnames, filenames in sorted(os.walk(model_results_dir)):
    for filename in filenames:
        if m := re.search(perm_results_regex, filename):
            model_name = m.group(1)
            print(model_name)
            _, cancer, analysis, target, data_type, model_code = model_name.split("_")

            figure_title = "{} {} ({})".format(
                cancer.upper(), target, model_code.upper()
            )

            perm_results = load(
                "{}/resp/{name}/{name}_perm_results.pkl".format(
                    model_results_dir, name=model_name
                )
            )

            colors = [
                (
                    "dark sky blue"
                    if data_type == "kraken"
                    else "burnt orange" if data_type == "star" else "purplish"
                )
            ]
            colors = sns.xkcd_palette(colors)

            fig, ax = plt.subplots(figsize=(fig_dim, fig_dim))
            perm_scores = perm_results["scores"]
            true_score = perm_results["true_score"]
            perm_pvalue = perm_results["pvalue"]
            # freedman-draconis rule
            bins = round(
                (np.max(perm_scores) - np.min(perm_scores))
                / (2 * iqr(perm_scores) / np.cbrt(perm_scores.size))
            )
            sns.histplot(
                perm_scores,
                bins=bins,
                kde=True,
                color=colors[0],
                stat="probability",
                edgecolor="white",
                lw=1,
                line_kws={"lw": 3},
            )
            ax.axvline(true_score, color="darkgrey", ls="--", lw=3)
            ax.set_title(
                figure_title, loc="center", pad=8, fontdict={"fontsize": title_fontsize}
            )
            ax.add_artist(
                AnchoredText(
                    # 'True AUROC = {:.2f}\n' r'$\itp$ = $\bf{:.{}}$'.format(
                    #     true_score, perm_pvalue, '2e' if perm_pvalue < 0.001
                    #     else '3f'),
                    "p = {:.{}}".format(
                        perm_pvalue, "2e" if perm_pvalue < 0.001 else "3f"
                    ),
                    loc="upper left",
                    frameon=False,
                    pad=0,
                    borderpad=0.2,
                    prop={"size": legend_fontsize},
                )
            )
            ax.set_xlabel("AUROC", fontsize=axis_fontsize)
            ax.set_ylabel("Probability", fontsize=axis_fontsize)
            ax.set_xticks(np.arange(0.0, 1.1, 0.2))
            ax.get_xaxis().set_major_formatter(
                ticker.FixedFormatter(["0", "0.2", "0.4", "0.6", "0.8", "1"])
            )
            ax.set_yticks(np.arange(0.0, 0.15, 0.02))
            ax.get_yaxis().set_major_formatter(
                ticker.FixedFormatter(
                    ["0", "0.02", "0.04", "0.06", "0.08", "0.10", "0.12", "0.14"]
                )
            )
            ax.set_xlim([0, 1])
            ax.set_ylim([0, 0.14])
            ax.tick_params(axis="both", labelsize=axis_fontsize)
            ax.tick_params(which="major", length=5, width=1.5)
            ax.tick_params(which="minor", width=1.5)
            plt.setp(ax.spines.values(), lw=1.5)
            ax.margins(0.01)
            ax.grid(False)
            ax.set_aspect(1.0 / ax.get_data_ratio())
            fig.tight_layout(pad=0.5, w_pad=0, h_pad=0)
            for fmt in args.file_format:
                fig.savefig(
                    "{}/{}_perm_hist.{}".format(args.out_dir, model_name, fmt),
                    format=fmt,
                    bbox_inches="tight",
                    # matplotlib GH#15497
                    dpi="figure" if fmt == "pdf" else fig_dpi,
                )
            pd.DataFrame(
                {
                    "perm_score": perm_scores,
                    "true_score": [true_score] + [""] * (len(perm_scores) - 1),
                    "p_value": [perm_pvalue] + [""] * (len(perm_scores) - 1),
                }
            ).to_csv(
                "{}/{}_perm_hist.tsv".format(args.out_dir, model_name),
                sep="\t",
                index=False,
            )
