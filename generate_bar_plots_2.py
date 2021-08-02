import os
import re
import sys
import warnings
from argparse import ArgumentParser

warnings.filterwarnings('ignore', category=FutureWarning,
                        module='rpy2.robjects.pandas2ri')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.api.types import is_object_dtype, is_string_dtype
import rpy2.rinterface_lib.embedded as r_embedded

r_embedded.set_initoptions(
    ('rpy2', '--quiet', '--no-save', '--max-ppsize=500000'))

import rpy2.robjects as robjects
import seaborn as sns
from joblib import dump, load
from matplotlib import ticker
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr

numpy2ri.activate()
pandas2ri.activate()

# suppress linux conda qt5 wayland warning
if sys.platform.startswith('linux'):
    os.environ['XDG_SESSION_TYPE'] = 'x11'

parser = ArgumentParser()
parser.add_argument('--results-dir', type=str, default='results',
                    help='results dir')
parser.add_argument('--out-dir', type=str, default='figures', help='out dir')
parser.add_argument('--filter', type=str, choices=['signif', 'all'],
                    default='signif',
                    help='response model filter')
parser.add_argument('--file-format', type=str, nargs='+',
                    choices=['png', 'pdf', 'svg', 'tif'], default=['png'],
                    help='save file format')
args = parser.parse_args()

model_results_dir = '{}/models'.format(args.results_dir)
analysis_results_dir = '{}/analysis'.format(args.results_dir)

os.makedirs(args.out_dir, mode=0o755, exist_ok=True)

data_types = ['kraken', 'htseq', 'combo']
metrics = ['roc_auc', 'pr_auc', 'balanced_accuracy']
metric_label = {'roc_auc': 'AUROC',
                'pr_auc': 'AUPRC',
                'balanced_accuracy': 'BCR'}
model_codes = ['rfe', 'lgr', 'edger', 'limma']
add_filters = pd.DataFrame(
    {'cancer': ['blca'], 'versus': ['cisplatin'], 'features': ['kraken']})

colors = ['crimson', 'ocean blue', 'greyish teal', 'steel grey']
colors = sns.xkcd_palette(colors)

title_fontsize = 16
x_axis_fontsize = 5 if args.filter == 'all' else 8
y_axis_fontsize = 12
legend_fontsize = 8
fig_let_fontsize = 48
fig_height = 4
fig_width = 10 if args.filter == 'all' else 6
fig_dpi = 300
bar_width = 0.8

plt.rcParams['figure.max_open_warning'] = 0
plt.rcParams['font.family'] = ['Nimbus Sans']

r_base = importr('base')

signif_hits = pd.read_csv(
    '{}/potential_hits.txt'.format(analysis_results_dir), sep='\t')
signif_hits = signif_hits.apply(
    lambda x: x.str.lower() if is_object_dtype(x) or is_string_dtype(x) else x)
signif_hits = signif_hits.loc[signif_hits['analysis'] == 'resp']
signif_hits = signif_hits.sort_values(
    by=['cancer', 'versus', 'features', 'how'])
signif_hits = signif_hits.loc[signif_hits.duplicated(
    subset=['cancer', 'versus', 'features'], keep=False)]

all_scores_dfs = {}
model_codes_regex = '|'.join(model_codes)
split_results_regex = re.compile(
    '^(.+?_(?:{}))_split_results\\.pkl$'.format(model_codes_regex))
for dirpath, dirnames, filenames in sorted(os.walk(model_results_dir)):
    for filename in filenames:
        if m := re.search(split_results_regex, filename):
            model_name = m.group(1)
            _, cancer, analysis, target, data_type, *rest = (
                model_name.split('_'))
            if data_type == 'htseq':
                model_code = '_'.join(rest[1:])
            else:
                model_code = '_'.join(rest)

            if (args.filter == 'signif'
                and not ((signif_hits['cancer'] == cancer)
                         & (signif_hits['versus'] == target)
                         & (signif_hits['features'] == data_type)).any()
                or ((add_filters['cancer'] == cancer)
                    & (add_filters['versus'] == target)
                    & (add_filters['features'] == data_type)).any()):
                continue

            split_results_file = '{}/{}'.format(dirpath, filename)
            print('Loading', split_results_file)
            split_results = load(split_results_file)
            for metric in metrics:
                scores = []
                for split_result in split_results:
                    if split_result is None:
                        scores.append(np.nan)
                    else:
                        scores.append(split_result['scores']['te'][metric])
                scores_df = pd.DataFrame({model_name: scores})
                if data_type not in all_scores_dfs:
                    all_scores_dfs[data_type] = {}
                if metric not in all_scores_dfs[data_type]:
                    all_scores_dfs[data_type][metric] = {}
                if model_code not in all_scores_dfs[data_type][metric]:
                    all_scores_dfs[data_type][metric][model_code] = scores_df
                else:
                    all_scores_dfs[data_type][metric][model_code] = pd.concat(
                        [all_scores_dfs[data_type][metric][model_code],
                         scores_df], axis=1)

for data_type in data_types:
    for metric in metrics:
        if data_type == 'kraken':
            data_type_label = 'Microbiome'
        elif data_type == 'htseq':
            data_type_label = 'Expression'
        else:
            data_type_label = 'Combo'
        x_pos = []
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=fig_dpi)
        for model_idx, model_code in enumerate(sorted(
                all_scores_dfs[data_type][metric].keys(),
                key=model_codes.index)):
            scores_df = all_scores_dfs[data_type][metric][model_code]
            x_pos.append(np.arange(0, scores_df.columns.size * 3, step=3)
                         if model_idx == 0 else
                         [x + bar_width for x in x_pos[-1]])
            mean_scores = scores_df.mean(axis=0)
            std_scores = scores_df.std(axis=0)
            upper_error = ((mean_scores + std_scores).clip(upper=1)
                           - mean_scores)
            lower_error = (mean_scores
                           - (mean_scores - std_scores).clip(lower=0))
            ax.bar(x_pos[-1], mean_scores, yerr=[lower_error, upper_error],
                   align='center', color=colors[model_idx], ecolor=colors[-1],
                   error_kw=dict(lw=0.75), label=model_code.upper(),
                   width=bar_width, zorder=3)
        ax.axhline(y=0.6, color=colors[-1], linestyle='--', lw=1, zorder=1)
        model_name_parts = pd.DataFrame(
            scores_df.columns. str.split('_', n=4).to_list(),
            columns=['program', 'cancer', 'analysis', 'target', 'rest'])
        x_labels = [' '.join(n) for n in zip(
            model_name_parts['cancer'].str.upper(),
            model_name_parts['target'].str.title())]
        ax.set_ylabel(metric_label[metric], fontsize=y_axis_fontsize,
                      labelpad=5)
        x_ticks = [x + bar_width for x in
                   np.arange(0, scores_df.columns.size * 3, step=3)]
        plt.xticks(x_ticks, x_labels, fontsize=x_axis_fontsize, rotation=60)
        ax.set_yticks(np.arange(0.0, 1.1, 0.2))
        ax.get_yaxis().set_major_formatter(ticker.FixedFormatter(
            ['0', '0.2', '0.4', '0.6', '0.8', '1']))
        ax.autoscale(axis='x', enable=None, tight=True)
        ax.set_ylim([0, 1.3])
        ax.tick_params(axis='x', direction='in', length=0, pad=1)
        ax.margins(0.01)
        ax.grid(False)
        legend = ax.legend(loc='upper right', labelspacing=0.25, frameon=False,
                           borderpad=0.1, handletextpad=0.25,
                           fontsize=legend_fontsize)
        legend.set_title(data_type_label, prop={'weight': 'bold',
                                                'size': y_axis_fontsize})
        legend._legend_box.align = 'right'
        renderer = fig.canvas.get_renderer()
        shift = max([text.get_window_extent(renderer).width
                     for text in legend.get_texts()])
        for text in legend.get_texts():
            text.set_ha('right')
            text.set_position((shift, 0))
        fig.tight_layout(pad=0.5, w_pad=0, h_pad=0)
        for fmt in args.file_format:
            fig.savefig('{}/{}_{}_bar_comp.{}'.format(args.out_dir, data_type,
                                                      metric, fmt),
                        format=fmt, bbox_inches='tight')
