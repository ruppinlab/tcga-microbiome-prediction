import os
import re
import sys
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pandas.api.types import is_object_dtype, is_string_dtype
import rpy2.rinterface_lib.embedded as r_embedded

r_embedded.set_initoptions(("rpy2", "--quiet", "--no-save", "--max-ppsize=500000"))

from joblib import load
from matplotlib import ticker
from matplotlib.container import BarContainer

# suppress linux conda qt5 wayland warning
if sys.platform.startswith("linux"):
    os.environ["XDG_SESSION_TYPE"] = "x11"
    os.environ["QT_QPA_PLATFORM"] = "xcb"

parser = ArgumentParser()
parser.add_argument("--results-dir", type=str, default="results", help="results dir")
parser.add_argument("--out-dir", type=str, default="figures/bar", help="out dir")
parser.add_argument(
    "--filter",
    type=str,
    choices=["signif", "all"],
    default="signif",
    help="response model filter",
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
analysis_results_dir = "{}/analysis".format(args.results_dir)

os.makedirs(args.out_dir, mode=0o755, exist_ok=True)

# for stripplot jitter
random_seed = 777
np.random.seed(random_seed)

data_types = ["kraken", "htseq", "combo"]
metrics = ["roc_auc", "pr_auc", "balanced_accuracy"]
metric_label = {"roc_auc": "AUROC", "pr_auc": "AUPRC", "balanced_accuracy": "BCR"}
model_codes = ["edge", "elgr", "srfe", "voom", "zinb"]

bcolors = ["#009E73", "#F0E442", "#0072B2"]
ecolor = "darkorange"

bcolors = sns.color_palette(bcolors)

title_fontsize = 16
x_axis_fontsize = 5 if args.filter == "all" else 10
y_axis_fontsize = 12
label_fontsize = 12
legend_fontsize = 12
fig_height = 4
fig_width = 10 if args.filter == "all" else 5
fig_dpi = 300
x_label_rotation = 60 if args.filter == "all" else 45

plt.rcParams["figure.max_open_warning"] = 0
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "Helvetica",
    "Nimbus Sans",
    "Arial",
    "DejaVu Sans",
    "sans-serif",
]

all_stats = pd.read_csv("{}/compared_runs.txt".format(analysis_results_dir), sep="\t")
all_stats = all_stats.apply(
    lambda x: x.str.lower() if is_object_dtype(x) or is_string_dtype(x) else x
)
all_stats = all_stats.loc[all_stats["analysis"] == "resp"]
all_stats = all_stats.sort_values(by=["cancer", "versus", "features", "how"])

signif_hits = pd.read_csv("{}/goodness_hits.txt".format(analysis_results_dir), sep="\t")
signif_hits = signif_hits.apply(
    lambda x: x.str.lower() if is_object_dtype(x) or is_string_dtype(x) else x
)
signif_hits = signif_hits.loc[signif_hits["analysis"] == "resp"]
signif_hits = signif_hits.sort_values(by=["cancer", "versus", "features", "how"])
signif_hits = signif_hits.loc[
    signif_hits.duplicated(subset=["cancer", "versus", "features"], keep=False)
]

score_dfs, p_adjs = {}, {}
model_codes_regex = "|".join(model_codes)
split_results_regex = re.compile(
    "^(.+?_(?:{}))_split_results\\.pkl$".format(model_codes_regex)
)
for dirpath, dirnames, filenames in sorted(os.walk(model_results_dir)):
    for filename in filenames:
        if m := re.search(split_results_regex, filename):
            model_name = m.group(1)
            _, cancer, analysis, target, data_type, *rest = model_name.split("_")
            if data_type == "htseq":
                model_code = "_".join(rest[1:])
            else:
                model_code = "_".join(rest)
            if (
                args.filter == "signif"
                and not (
                    (signif_hits["cancer"] == cancer)
                    & (signif_hits["versus"] == target)
                    & (signif_hits["features"] == data_type)
                ).any()
            ):
                continue
            if data_type not in p_adjs:
                p_adjs[data_type] = {}
            if model_code not in p_adjs[data_type]:
                p_adjs[data_type][model_code] = []
            p_adj = all_stats.loc[
                (all_stats["cancer"] == cancer)
                & (all_stats["versus"] == target)
                & (all_stats["features"] == data_type)
                & (all_stats["how"] == model_code),
                "p_adj",
            ].item()
            p_adjs[data_type][model_code].append(p_adj)
            split_results_file = "{}/{}".format(dirpath, filename)
            print("Loading", split_results_file)
            split_results = load(split_results_file)
            for metric in metrics:
                scores = []
                for split_result in split_results:
                    if split_result is None:
                        scores.append(np.nan)
                    else:
                        scores.append(split_result["scores"]["te"][metric])
                if data_type not in score_dfs:
                    score_dfs[data_type] = {}
                score_df = pd.DataFrame(
                    {
                        "cancer": cancer,
                        "target": target,
                        "type": model_code,
                        "score": scores,
                    }
                )
                if metric not in score_dfs[data_type]:
                    score_dfs[data_type][metric] = score_df
                else:
                    score_dfs[data_type][metric] = pd.concat(
                        [score_dfs[data_type][metric], score_df], axis=0
                    )

for data_type in data_types:
    for metric in metrics:
        score_df = score_dfs[data_type][metric].copy()
        score_df["model"] = pd.Categorical(
            score_df["cancer"].str.upper() + " " + score_df["target"], ordered=True
        )
        score_df["type"] = pd.Categorical(
            score_df["type"].str.upper(),
            categories=[m.upper() for m in model_codes],
            ordered=True,
        )
        score_df["type"] = score_df["type"].cat.remove_unused_categories()
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=fig_dpi)
        sns.barplot(
            x="model",
            y="score",
            hue="type",
            data=score_df,
            ci="sd",
            capsize=0.1,
            palette=bcolors,
            errcolor=ecolor,
            errwidth=1.75,
            saturation=1,
        )
        # bar_labels = ['***' if p <= 0.0001 else '**' if p <= 0.001 else
        #               '*' if p <= 0.01 else '$^{ns}$' for p in np.ravel([
        #                   p_adjs[data_type][k] for k in
        #                   score_df['type'].cat.categories.str.lower()])]
        # for bar, label in zip(ax.patches, bar_labels):
        #     ax.annotate(label, (bar.get_x() + bar.get_width() / 2, 1.0),
        #                 ha='center', va='bottom', size=label_fontsize,
        #                 xytext=(0, 1), textcoords='offset points')
        for line in ax.get_lines():
            x, y = line.get_data()
            line.set_data(x, np.clip(y, 0, 1))
        sns.stripplot(
            x="model",
            y="score",
            hue="type",
            data=score_df,
            palette=bcolors,
            edgecolor="black",
            dodge=True,
            jitter=0.15,
            size=2.5,
            linewidth=0.8,
            alpha=0.6,
            zorder=2,
        )
        # ax.axhline(y=0.6, color='darkgrey', linestyle='--', lw=2, zorder=0)
        ax.autoscale(axis="x", enable=None, tight=True)
        ax.tick_params(which="major", length=3, width=1.25)
        ax.tick_params(
            axis="x",
            direction="out",
            labelsize=x_axis_fontsize,
            labelrotation=x_label_rotation,
            length=3,
            width=1.25,
            pad=0,
        )
        ax.set_xlabel(None)
        ax.set_ylabel(metric_label[metric], fontsize=y_axis_fontsize, labelpad=5)
        ax.set_yticks(np.arange(0.0, 1.1, 0.2))
        ax.get_yaxis().set_major_formatter(
            ticker.FixedFormatter(["0", "0.2", "0.4", "0.6", "0.8", "1"])
        )
        ax.set_ylim([-0.01, 1.15])
        plt.setp(ax.spines.values(), lw=1.25)
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        ax.spines.left.set_bounds(-0.01, 1)
        ax.margins(0.01)
        ax.grid(True, alpha=0.3, axis="y", which="major")
        ax.set_axisbelow(True)
        handles, labels = ax.get_legend_handles_labels()
        handles, labels = zip(
            *((h, l) for h, l in zip(handles, labels) if isinstance(h, BarContainer))
        )
        legend = ax.legend(
            handles=handles,
            labels=labels,
            loc="upper right",
            labelspacing=0.25,
            frameon=False,
            borderpad=0,
            handletextpad=0.25,
            fontsize=legend_fontsize,
            ncol=len(ax.lines),
            columnspacing=1,
        )
        # legend.set_title('Microbiome' if data_type == 'kraken' else
        #                  'Expression' if data_type == 'htseq' else
        #                  'Combo', prop={'weight': 'bold',
        #                                 'size': y_axis_fontsize})
        legend._legend_box.align = "right"
        fig.tight_layout(pad=0.5, w_pad=0, h_pad=0)
        for fmt in args.file_format:
            fig.savefig(
                "{}/{}_{}_bar_comp.{}".format(args.out_dir, data_type, metric, fmt),
                format=fmt,
                bbox_inches="tight",
            )
        score_df.to_csv(
            "{}/{}_{}_bar_comp.tsv".format(args.out_dir, data_type, metric),
            sep="\t",
            index=False,
        )
