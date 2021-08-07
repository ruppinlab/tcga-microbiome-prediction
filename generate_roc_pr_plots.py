import os
import re
import sys
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from joblib import load
from matplotlib import ticker

# suppress linux conda qt5 wayland warning
if sys.platform.startswith('linux'):
    os.environ['XDG_SESSION_TYPE'] = 'x11'

parser = ArgumentParser()
parser.add_argument('--results-dir', type=str, default='results',
                    help='results dir')
parser.add_argument('--out-dir', type=str, default='figures', help='out dir')
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
legend_fontsize = 10
fig_let_fontsize = 48
fig_dim = 4
fig_dpi = 300

plt.rcParams['figure.max_open_warning'] = 0
plt.rcParams['font.family'] = [
    'Nimbus Sans' if sys.platform.startswith('linux') else 'sans']

model_codes_regex = '|'.join(args.model_code)
split_results_regex = re.compile(
    '^(.+?_(?:{}))_split_results\\.pkl$'.format(model_codes_regex))
for dirpath, dirnames, filenames in sorted(os.walk(model_results_dir)):
    for filename in filenames:
        if m := re.search(split_results_regex, filename):
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

            split_results = []
            split_results.append(load(
                '{}/resp/{name}/{name}_split_results.pkl'
                .format(model_results_dir, name=model_name)))
            if data_type in ('kraken', 'htseq'):
                dataset_name = '_'.join(model_name.split('_')[:-1])
                clinical_model_name = '_'.join(
                    [dataset_name, 'svm' if model_code in ('rfe') else 'lgr',
                     'clinical'])
                split_results.append(
                    load('{}/resp/{name}/{name}_split_results.pkl'
                         .format(model_results_dir, name=clinical_model_name)))
            else:
                for new_data_type in ('htseq_counts', 'kraken'):
                    new_model_code = (
                        'edger' if new_data_type == 'htseq_counts'
                        and model_code == 'limma' else rest[-1])
                    new_model_name = '_'.join(
                        model_name.split('_')[:-2]
                        + [new_data_type, new_model_code])
                    split_results.append(load(
                        '{}/resp/{name}/{name}_split_results.pkl'
                        .format(model_results_dir, name=new_model_name)))

            if data_type == 'kraken':
                colors = ['dark sky blue', 'purplish']
            elif data_type == 'htseq':
                colors = ['burnt orange', 'turquoise']
            else:
                colors = ['purplish', 'burnt orange', 'dark sky blue']

            colors.append('steel grey')
            colors = sns.xkcd_palette(colors)

            # roc curves
            fig, ax = plt.subplots(figsize=(fig_dim, fig_dim), dpi=fig_dpi)
            for ridx, _ in enumerate(split_results):
                tprs, roc_scores = [], []
                mean_fpr = np.linspace(0, 1, 1000)
                for split_result in split_results[ridx]:
                    if split_result is None:
                        continue
                    tprs.append(np.interp(
                        mean_fpr, split_result['scores']['te']['fpr'],
                        split_result['scores']['te']['tpr']))
                    tprs[-1][0] = 0.0
                    roc_scores.append(split_result['scores']['te']['roc_auc'])
                mean_tpr = np.mean(tprs, axis=0)
                mean_tpr[-1] = 1.0
                std_tpr = np.std(tprs, axis=0)
                tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
                tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
                if data_type == 'combo':
                    dtype_label = ('Combo' if ridx == 0 else
                                   'Expression' if ridx == 1 else 'Microbiome')
                    label = '{} + Clinical'.format(dtype_label)
                    color = colors[ridx]
                    zorder = 2.5 if ridx == 0 else 2.2 if ridx == 1 else 2
                else:
                    if ridx == 0:
                        label = '{} + Clinical'.format(data_type_label)
                        color = colors[0]
                        zorder = 2.5
                    else:
                        label = 'Clinical'
                        color = colors[-1]
                        zorder = 2
                ax.plot(mean_fpr, mean_tpr, alpha=0.8, color=color, lw=2,
                        label=r'{} AUROC = $\bf{{{:.2f}}}$'.format(
                            label, np.mean(roc_scores)), zorder=zorder)
                ax.fill_between(mean_fpr, tprs_lower, tprs_upper, alpha=0.1,
                                color=color, zorder=zorder)
            ax.plot([0, 1], [0, 1], alpha=0.2, color='darkgrey',
                    linestyle='--', lw=1.5, zorder=1)
            ax.set_xlabel('False positive rate', fontsize=axis_fontsize,
                          labelpad=5)
            ax.set_ylabel('True positive rate', fontsize=axis_fontsize,
                          labelpad=5)
            ax.set_yticks(np.arange(0.0, 1.1, 0.2))
            ax.get_yaxis().set_major_formatter(ticker.FixedFormatter(
                ['0', '0.2', '0.4', '0.6', '0.8', '1']))
            ax.set_xticks(np.arange(0.0, 1.1, 0.2))
            ax.get_xaxis().set_major_formatter(ticker.FixedFormatter(
                ['0', '0.2', '0.4', '0.6', '0.8', '1']))
            ax.set_xlim([-0.01, 1.01])
            ax.set_ylim([-0.01, 1.01])
            ax.tick_params(axis='both', labelsize=axis_fontsize)
            ax.tick_params(which='major', width=1)
            ax.tick_params(which='major', length=5)
            ax.tick_params(which='minor', width=1)
            ax.margins(0)
            ax.grid(False)
            legend = ax.legend(loc='lower right', frameon=False, borderpad=0.1,
                               prop={'size': legend_fontsize})
            legend.set_title(legend_title, prop={'weight': 'bold',
                                                 'size': axis_fontsize})
            legend._legend_box.align = 'right'
            for item in legend.legendHandles:
                item.set_visible(False)
            renderer = fig.canvas.get_renderer()
            text_widths = [text.get_window_extent(renderer).width
                           for text in legend.get_texts()]
            max_width = max(text_widths)
            shifts = [max_width - w for w in text_widths]
            for i, text in enumerate(legend.get_texts()):
                text.set_ha('right')
                text.set_position((shifts[i], 0))
            ax.set_aspect(1.0 / ax.get_data_ratio())
            fig.tight_layout(pad=0.5, w_pad=0, h_pad=0)
            for fmt in args.file_format:
                fig.savefig('{}/{}_roc_auc.{}'.format(args.out_dir, model_name,
                                                      fmt),
                            format=fmt, bbox_inches='tight')

            # pr curves
            fig, ax = plt.subplots(figsize=(fig_dim, fig_dim), dpi=fig_dpi)
            for ridx, _ in enumerate(split_results):
                pres, pr_scores = [], []
                mean_rec = np.linspace(0, 1, 1000)
                for split_result in split_results[ridx]:
                    if split_result is None:
                        continue
                    pres.append(np.interp(
                        mean_rec, split_result['scores']['te']['rec'][::-1],
                        split_result['scores']['te']['pre'][::-1]))
                    pr_scores.append(split_result['scores']['te']['pr_auc'])
                mean_pre = np.mean(pres, axis=0)
                std_pre = np.std(pres, axis=0)
                pres_upper = np.minimum(mean_pre + std_pre, 1)
                pres_lower = np.maximum(mean_pre - std_pre, 0)
                if data_type == 'combo':
                    dtype_label = ('Combo' if ridx == 0 else
                                   'Expression' if ridx == 1 else 'Microbiome')
                    label = '{} + Clinical'.format(dtype_label)
                    color = colors[ridx]
                    zorder = 2.5 if ridx == 0 else 2.2 if ridx == 1 else 2
                else:
                    if ridx == 0:
                        label = '{} + Clinical'.format(data_type_label)
                        color = colors[1]
                        zorder = 2.5
                    else:
                        label = 'Clinical'
                        color = colors[-1]
                        zorder = 2
                ax.step(mean_rec, mean_pre, alpha=0.8, color=color, lw=2,
                        label=r'{} AUPRC = $\bf{{{:.2f}}}$'.format(
                            label, np.mean(pr_scores)), where='post',
                        zorder=zorder)
                ax.fill_between(mean_rec, pres_lower, pres_upper, alpha=0.1,
                                color=color, zorder=zorder)
            ax.set_xlabel('Recall', fontsize=axis_fontsize, labelpad=5)
            ax.set_ylabel('Precision', fontsize=axis_fontsize, labelpad=5)
            ax.set_xticks(np.arange(0.0, 1.1, 0.2))
            ax.get_xaxis().set_major_formatter(ticker.FixedFormatter(
                ['0', '0.2', '0.4', '0.6', '0.8', '1']))
            ax.set_yticks(np.arange(0.0, 1.1, 0.2))
            ax.get_yaxis().set_major_formatter(ticker.FixedFormatter(
                ['0', '0.2', '0.4', '0.6', '0.8', '1']))
            ax.set_xlim([-0.01, 1.01])
            ax.set_ylim([-0.01, 1.01])
            ax.tick_params(axis='both', labelsize=axis_fontsize)
            ax.tick_params(which='major', width=1)
            ax.tick_params(which='major', length=5)
            ax.tick_params(which='minor', width=1)
            ax.margins(0)
            ax.grid(False)
            legend = ax.legend(loc='lower right', frameon=False, borderpad=0.1,
                               prop={'size': legend_fontsize})
            legend.set_title(legend_title, prop={'weight': 'bold',
                                                 'size': axis_fontsize})
            legend._legend_box.align = 'right'
            for item in legend.legendHandles:
                item.set_visible(False)
            renderer = fig.canvas.get_renderer()
            text_widths = [text.get_window_extent(renderer).width
                           for text in legend.get_texts()]
            max_width = max(text_widths)
            shifts = [max_width - w for w in text_widths]
            for i, text in enumerate(legend.get_texts()):
                text.set_ha('right')
                text.set_position((shifts[i], 0))
            ax.set_aspect(1.0 / ax.get_data_ratio())
            fig.tight_layout(pad=0.5, w_pad=0, h_pad=0)
            for fmt in args.file_format:
                fig.savefig('{}/{}_pr_auc.{}'.format(args.out_dir, model_name,
                                                     fmt),
                            format=fmt, bbox_inches='tight')
