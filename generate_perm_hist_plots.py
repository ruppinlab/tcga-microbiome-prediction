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

title_fontsize = 16
axis_fontsize = 12
legend_fontsize = 9
fig_let_fontsize = 48
fig_dim = 4
fig_dpi = 300

plt.rcParams['figure.max_open_warning'] = 0
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Nimbus Sans', 'DejaVu Sans', 'sans']

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
            legend_title = '{} {} ({})'.format(cancer.upper(), target.title(),
                                               model_code.upper())
            data_type_label = ('Expression' if data_type == 'htseq' else
                               'Microbiome')

            perm_results = load('{}/resp/{name}/{name}_perm_results.pkl'
                                .format(model_results_dir, name=model_name))

            color = sns.xkcd_palette([
                'dark sky blue' if data_type == 'kraken'
                else 'burnt_orange' if data_type == 'htseq'
                else 'purplish'])[0]

            fig, ax = plt.subplots(figsize=(fig_dim, fig_dim), dpi=fig_dpi)
            perm_scores = perm_results['scores']
            true_score = perm_results['true_score']
            perm_pvalue = perm_results['pvalue']
            # freedman-draconis rule
            bins = round((np.max(perm_scores) - np.min(perm_scores))
                         / (2 * iqr(perm_scores) / np.cbrt(perm_scores.size)))
            sns.histplot(perm_scores, bins=bins, kde=True, color=color,
                         stat='probability', edgecolor='white')
            ax.axvline(true_score, ls='--', color='darkgrey')
            ax.add_artist(AnchoredText(
                 r'$\bf{{{}}}$' '\n' r'True AUROC = {:.2f}' '\n'
                 r'$\itp$ = {:.{}}'.format(
                     legend_title.replace(' ', '\\ '), true_score,
                     perm_pvalue, '2e' if perm_pvalue < 0.001 else '3f'),
                 loc='upper left', frameon=False, pad=0,
                 prop={'size': legend_fontsize}))
            ax.set_xlabel('AUROC', fontsize=axis_fontsize)
            ax.set_ylabel('Probability', fontsize=axis_fontsize)
            ax.set_xticks(np.arange(0.0, 1.1, 0.2))
            ax.get_xaxis().set_major_formatter(ticker.FixedFormatter(
                ['0', '0.2', '0.4', '0.6', '0.8', '1']))
            ax.set_yticks(np.arange(0.0, 0.3, 0.05))
            ax.get_yaxis().set_major_formatter(ticker.FixedFormatter(
                ['0', '0.05', '0.1', '0.15', '0.2', '0.25']))
            ax.set_xlim([0, 1])
            ax.set_ylim([0, 0.275])
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
