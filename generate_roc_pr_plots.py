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
parser.add_argument('--out-dir', type=str, default='figures/roc_pr',
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
legend_fontsize = 10
fig_let_fontsize = 48
fig_dim = 4
fig_dpi = 300

plt.rcParams['figure.max_open_warning'] = 0
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Nimbus Sans', 'DejaVu Sans', 'sans']

pipe_step_type_regex = re.compile(
    r'^({})\d+$'.format('|'.join(['slr', 'trf', 'clf'])))

param_types = {'edger': ['slr__k'],
               'lgr': ['slr__estimator__C', 'slr__estimator__l1_ratio'],
               'limma': ['slr__k'],
               'rfe': ['clf__n_features_to_select']}

metrics = ['roc_auc', 'average_precision']
metric_labels = ['AUROC', 'AVPRE']

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
                mean_fprs = np.linspace(0, 1, 1000)
                for split_result in split_results[ridx]:
                    if split_result is None:
                        continue
                    tprs.append(np.interp(
                        mean_fprs, split_result['scores']['te']['fpr'],
                        split_result['scores']['te']['tpr']))
                    tprs[-1][0] = 0.0
                    roc_scores.append(split_result['scores']['te']['roc_auc'])
                mean_tprs = np.mean(tprs, axis=0)
                mean_tprs[-1] = 1.0
                std_tprs = np.std(tprs, axis=0)
                tprs_upper = np.minimum(mean_tprs + std_tprs, 1)
                tprs_lower = np.maximum(mean_tprs - std_tprs, 0)
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
                ax.plot(mean_fprs, mean_tprs, alpha=0.8, color=color, lw=2,
                        label=r'{} AUROC = $\bf{:.2f}$'.format(
                            label, np.mean(roc_scores)), zorder=zorder)
                ax.fill_between(mean_fprs, tprs_lower, tprs_upper, alpha=0.1,
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
                mean_recs = np.linspace(0, 1, 1000)
                for split_result in split_results[ridx]:
                    if split_result is None:
                        continue
                    pres.append(np.interp(
                        mean_recs, split_result['scores']['te']['rec'][::-1],
                        split_result['scores']['te']['pre'][::-1]))
                    pr_scores.append(split_result['scores']['te']['pr_auc'])
                mean_pres = np.mean(pres, axis=0)
                std_pres = np.std(pres, axis=0)
                pres_upper = np.minimum(mean_pres + std_pres, 1)
                pres_lower = np.maximum(mean_pres - std_pres, 0)
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
                ax.step(mean_recs, mean_pres, alpha=0.8, color=color, lw=2,
                        label=r'{} AUPRC = $\bf{:.2f}$'.format(
                            label, np.mean(pr_scores)), where='post',
                        zorder=zorder)
                ax.fill_between(mean_recs, pres_lower, pres_upper, alpha=0.1,
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

            # num selected features vs scores
            param_cv_scores = load('{}/resp/{name}/{name}_param_cv_scores.pkl'
                                   .format(model_results_dir, name=model_name))

            if data_type == 'combo':
                colors = ['indigo', 'magenta']
                colors = sns.xkcd_palette(colors)

            for param in param_cv_scores:
                param_parts = param.split('__')
                param_parts_start_idx = [i for i, p in enumerate(param_parts)
                                         if pipe_step_type_regex.match(p)][-1]
                param_parts[param_parts_start_idx] = pipe_step_type_regex.sub(
                    r'\1', param_parts[param_parts_start_idx])
                param_type = '__'.join(param_parts[param_parts_start_idx:])
                if param_type not in param_types[model_code]:
                    continue
                fig, ax = plt.subplots(figsize=(fig_dim, fig_dim), dpi=fig_dpi)
                mean_cv_scores, std_cv_scores = {}, {}
                for metric in metrics:
                    param_metric_scores = (
                        param_cv_scores[param][metric]['scores'])
                    param_metric_stdev = (
                        param_cv_scores[param][metric]['stdev'])
                    if any(len(scores) > 1 for scores in param_metric_scores):
                        mean_cv_scores[metric], std_cv_scores[metric] = [], []
                        for param_value_scores in param_metric_scores:
                            mean_cv_scores[metric].append(
                                np.mean(param_value_scores))
                            std_cv_scores[metric].append(
                                np.std(param_value_scores))
                    else:
                        mean_cv_scores[metric] = np.ravel(param_metric_scores)
                        std_cv_scores[metric] = np.ravel(param_metric_stdev)
                if model_code in ('edger', 'limma', 'rfe'):
                    x_axis = np.insert(np.linspace(2, 400, num=200, dtype=int),
                                       0, 1)
                    x_label = 'Number of selected features'
                    ax.set_xlim([0, max(x_axis)])
                    ax.set_xticks([1] + list(range(50, 450, 50)))
                    param_ext = 'k'
                elif param_parts[-1] == 'C':
                    x_axis = (np.logspace(-2, 3, 6) if data_type == 'kraken'
                              else np.logspace(-2, 1, 4))
                    x_label = 'C'
                    ax.set_xlim([min(x_axis), max(x_axis)])
                    ax.set_xscale('log')
                    ax.set_xticks(x_axis)
                    param_ext = 'c'
                elif param_parts[-1] == 'l1_ratio':
                    x_axis = np.array([0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95,
                                       0.99, 1.])
                    x_label = 'L1 ratio'
                    ax.set_xlim([min(x_axis), max(x_axis)])
                    ax.set_xticks(x_axis)
                    ax.get_xaxis().set_major_formatter(ticker.FixedFormatter(
                        ['0.1', '0.3', '0.5', '0.7', '0.8', '0.9', '', '',
                         '1']))
                    param_ext = 'l1r'
                for metric_idx, metric in enumerate(metrics):
                    zorder = (2.5 if metric_idx == 0 else
                              2.2 if metric_idx == 1 else 2)
                    ax.plot(x_axis, mean_cv_scores[metric],
                            color=colors[metric_idx], lw=2, alpha=0.8,
                            label='{}'.format(metric_labels[metric_idx]),
                            zorder=zorder)
                    ax.fill_between(
                        x_axis,
                        [m - s for m, s in zip(mean_cv_scores[metric],
                                               std_cv_scores[metric])],
                        [m + s for m, s in zip(mean_cv_scores[metric],
                                               std_cv_scores[metric])],
                        alpha=0.1, color=colors[metric_idx], zorder=zorder)
                ax.set_xlabel(x_label, fontsize=axis_fontsize)
                ax.set_ylabel('Score', fontsize=axis_fontsize)
                ax.set_ylim([0.0, 1.0])
                ax.set_yticks(np.arange(0.0, 1.1, 0.2))
                ax.get_yaxis().set_major_formatter(ticker.FixedFormatter(
                    ['0', '0.2', '0.4', '0.6', '0.8', '1']))
                ax.tick_params(axis='both', labelsize=axis_fontsize)
                ax.tick_params(which='major', width=1)
                ax.tick_params(which='major', length=5)
                ax.tick_params(which='minor', width=1)
                ax.margins(0)
                ax.grid(True, alpha=0.3)
                legend = ax.legend(loc='lower right', borderpad=0.2,
                                   prop={'size': legend_fontsize})
                legend.set_title(legend_title, prop={'weight': 'bold',
                                                     'size': axis_fontsize})
                legend._legend_box.align = 'right'
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
                    fig.savefig('{}/{}_{}_vs_score.{}'.format(
                        args.out_dir, model_name, param_ext, fmt),
                                format=fmt, bbox_inches='tight')
