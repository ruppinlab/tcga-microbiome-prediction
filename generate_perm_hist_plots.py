import os
import re
import sys
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from joblib import load
from matplotlib import ticker
from matplotlib.offsetbox import AnchoredText
from scipy.stats import iqr

# suppress linux conda qt5 wayland warning
if sys.platform.startswith('linux'):
    os.environ['XDG_SESSION_TYPE'] = 'x11'

parser = ArgumentParser()
parser.add_argument('--results-dir', type=str, default='results',
                    help='results dir')
parser.add_argument('--out-dir', type=str, default='figures/perm_hist',
                    help='out dir')
parser.add_argument('--model-code', type=str, nargs='+',
                    choices=['edger', 'lgr', 'limma', 'rfe'],
                    default=['edger', 'lgr', 'limma', 'rfe'],
                    help='response model code filter')
parser.add_argument('--file-format', type=str, nargs='+',
                    choices=['png', 'pdf', 'svg', 'tif'], default=['png'],
                    help='save file format')
args = parser.parse_args()

model_results_dir = '{}/models'.format(args.results_dir)

os.makedirs(args.out_dir, mode=0o755, exist_ok=True)

title_fontsize = 14
axis_fontsize = 12
legend_fontsize = 12
fig_let_fontsize = 48
fig_dim = 4
fig_dpi = 300

plt.rcParams['figure.max_open_warning'] = 0
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'sans']

model_codes_regex = '|'.join(args.model_code)
perm_results_regex = re.compile(
    '^(.+?_(?:{}))_perm_results\\.pkl$'.format(model_codes_regex))
for dirpath, dirnames, filenames in sorted(os.walk(model_results_dir)):
    for filename in filenames:
        if m := re.search(perm_results_regex, filename):
            model_name = m.group(1)
            print(model_name)
            _, cancer, analysis, target, data_type, *rest = (
                model_name.split('_'))
            if data_type == 'htseq':
                model_code = '_'.join(rest[1:])
            else:
                model_code = '_'.join(rest)
            figure_title = '{} {} ({})'.format(cancer.upper(), target,
                                               model_code.upper())
            data_type_label = ('Expression' if data_type == 'htseq' else
                               'Microbiome')

            perm_results = load('{}/resp/{name}/{name}_perm_results.pkl'
                                .format(model_results_dir, name=model_name))

            colors = ['dark sky blue' if data_type == 'kraken' else
                      'burnt orange' if data_type == 'htseq' else 'purplish']
            colors = sns.xkcd_palette(colors)

            fig, ax = plt.subplots(figsize=(fig_dim, fig_dim), dpi=fig_dpi)
            perm_scores = perm_results['scores']
            true_score = perm_results['true_score']
            perm_pvalue = perm_results['pvalue']
            # freedman-draconis rule
            bins = round((np.max(perm_scores) - np.min(perm_scores))
                         / (2 * iqr(perm_scores) / np.cbrt(perm_scores.size)))
            sns.histplot(perm_scores, bins=bins, kde=True, color=colors[0],
                         stat='probability', edgecolor='white')
            ax.axvline(true_score, ls='--', color='darkgrey')
            ax.set_title(figure_title, loc='left', y=1.0, pad=4,
                         fontdict={'fontsize': title_fontsize})
            ax.add_artist(AnchoredText(
                # r'True AUROC = {:.2f}' '\n' r'$\itp$ = $\bf{:.{}}$'.format(
                #     true_score, perm_pvalue, '2e' if perm_pvalue < 0.001
                #     else '3f'),
                r'p = {:.{}}'.format(
                    perm_pvalue, '2e' if perm_pvalue < 0.001 else '3f'),
                loc='upper left', frameon=False, pad=0,
                prop={'size': legend_fontsize}))
            ax.set_xlabel('AUROC', fontsize=axis_fontsize)
            ax.set_ylabel('Probability', fontsize=axis_fontsize)
            ax.set_xticks(np.arange(0.0, 1.1, 0.2))
            ax.get_xaxis().set_major_formatter(ticker.FixedFormatter(
                ['0', '0.2', '0.4', '0.6', '0.8', '1']))
            ax.set_yticks(np.arange(0.0, 0.15, 0.02))
            ax.get_yaxis().set_major_formatter(ticker.FixedFormatter(
                ['0', '0.02', '0.04', '0.06', '0.08', '0.10', '0.12', '0.14']))
            ax.set_xlim([0, 1])
            ax.set_ylim([0, 0.14])
            ax.tick_params(axis='both', labelsize=axis_fontsize)
            ax.tick_params(which='major', width=1)
            ax.tick_params(which='major', length=5)
            ax.tick_params(which='minor', width=1)
            # ax.margins(0)
            ax.grid(False)
            ax.set_aspect(1.0 / ax.get_data_ratio())
            fig.tight_layout(pad=0.5, w_pad=0, h_pad=0)
            for fmt in args.file_format:
                fig.savefig('{}/{}_perm_hist.{}'.format(args.out_dir,
                                                        model_name, fmt),
                            format=fmt, bbox_inches='tight')
            tsv_filename = '{}/{}_perm_hist.tsv'.format(
                args.out_dir, model_name
            )
            with open(tsv_filename, "w") as fh:
                print(
                    "cancer", "target", "data_type", "model_code",
                    "replicate", "score", sep="\t", file=fh,
                )
                for i, score in enumerate(perm_scores):
                    print(
                        cancer, target, data_type_label,
                        model_code.upper(), i + 1, score,
                        sep="\t", file=fh,
                    )
                print(
                    cancer, target, data_type_label,
                    model_code.upper(), "True_score", true_score,
                    sep="\t", file=fh,
                )
