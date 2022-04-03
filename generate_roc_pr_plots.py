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

title_fontsize = 20
axis_fontsize = 18
legend_fontsize = 18
fig_dim = 4
fig_dpi = 300

plt.rcParams['figure.max_open_warning'] = 0
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Nimbus Sans', 'Arial',
                                   'DejaVu Sans', 'sans-serif']

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

            dtype_labels = []
            dtype_labels.append('Combo' if data_type == 'combo' else
                                'Expression' if data_type == 'htseq' else
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
                for new_data_type in ('htseq', 'kraken'):
                    dtype_labels.append((
                        'Expression' if new_data_type == 'htseq' else
                        'Microbiome'))
                    new_model_code = ('edger' if new_data_type == 'htseq'
                                      and model_code == 'limma' else rest[-1])
                    new_model_name_parts = model_name.split('_')[:-2]
                    new_model_name_parts.append(new_data_type)
                    if new_data_type == 'htseq':
                        new_model_name_parts.append('counts')
                    new_model_name_parts.append(new_model_code)
                    new_model_name = '_'.join(new_model_name_parts)
                    split_results.append(load(
                        '{}/resp/{name}/{name}_split_results.pkl'
                        .format(model_results_dir, name=new_model_name)))

            dtype_labels.append('Clinical')
            abbr_dtype_labels = ['Express' if l == 'Expression' else
                                 'Microbe' if l == 'Microbiome' else
                                 l for l in dtype_labels]

            figure_title = '{} {} ({})'.format(cancer.upper(), target,
                                               model_code.upper())

            # roc curves
            if data_type == 'kraken':
                colors = ['dark sky blue', 'steel grey']
            elif data_type == 'htseq':
                colors = ['burnt orange', 'steel grey']
            else:
                colors = ['purplish', 'burnt orange', 'dark sky blue']

            colors = sns.xkcd_palette(colors)

            tsv_scores = {k: [] for k in ['dtype', 'split', 'fpr', 'tpr']}
            fig, ax = plt.subplots(figsize=(fig_dim, fig_dim))
            for ridx, _ in enumerate(split_results):
                tprs, roc_scores = [], []
                mean_fprs = np.linspace(0, 1, 1000)
                for split_idx, split_result in enumerate(split_results[ridx]):
                    if split_result is None:
                        continue
                    fpr = split_result['scores']['te']['fpr']
                    tpr = split_result['scores']['te']['tpr']
                    tprs.append(np.interp(mean_fprs, fpr, tpr))
                    tprs[-1][0] = 0.0
                    roc_scores.append(split_result['scores']['te']['roc_auc'])
                    if ridx == 0:
                        tsv_data_type = data_type
                    elif data_type == 'combo':
                        tsv_data_type = 'htseq' if ridx == 1 else 'kraken'
                    else:
                        tsv_data_type = 'clinical'
                    tsv_scores['dtype'].extend([tsv_data_type] * len(fpr))
                    tsv_scores['split'].extend([split_idx + 1] * len(fpr))
                    tsv_scores['fpr'].extend(fpr)
                    tsv_scores['tpr'].extend(tpr)
                mean_tprs = np.mean(tprs, axis=0)
                mean_tprs[-1] = 1.0
                std_tprs = np.std(tprs, axis=0)
                tprs_upper = np.minimum(mean_tprs + std_tprs, 1)
                tprs_lower = np.maximum(mean_tprs - std_tprs, 0)
                if data_type == 'combo':
                    label = '+'.join([dtype_labels[ridx], dtype_labels[-1]])
                    zorder = 2.5 if ridx == 0 else 2.2 if ridx == 1 else 2
                elif ridx == 0:
                    label = '+'.join([dtype_labels[ridx], dtype_labels[-1]])
                    zorder = 2.5
                else:
                    label = dtype_labels[-1]
                    zorder = 2
                ax.plot(mean_fprs, mean_tprs, alpha=0.8, color=colors[ridx],
                        label=('AUROC = {:.2f}'.format(np.mean(roc_scores))
                               if ridx == 0 else None), lw=2, zorder=zorder)
                ax.fill_between(mean_fprs, tprs_lower, tprs_upper, alpha=0.1,
                                color=colors[ridx], zorder=zorder)
            ax.plot([0, 1], [0, 1], alpha=0.2, color='darkgrey',
                    linestyle='--', lw=1.5, zorder=1)
            ax.set_title(figure_title, loc='center', pad=8,
                         fontdict={'fontsize': title_fontsize})
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
            ax.tick_params(which='major', length=5, width=1)
            ax.tick_params(which='minor', width=1)
            ax.margins(0)
            ax.grid(False)
            legend = ax.legend(loc='lower right', borderpad=0.1,
                               borderaxespad=0.1, frameon=False,
                               labelspacing=0.2, fontsize=legend_fontsize)
            legend.set_title(
                '+'.join([abbr_dtype_labels[0], abbr_dtype_labels[-1]]),
                prop={'weight': 'regular', 'size': legend_fontsize})
            legend._legend_box.align = 'right'
            for item in legend.legendHandles:
                item.set_visible(False)
            text_widths = [
                text.get_window_extent(fig.canvas.get_renderer()).width
                for text in legend.get_texts()]
            max_width = max(text_widths)
            shifts = [max_width - w for w in text_widths]
            for i, text in enumerate(legend.get_texts()):
                text.set_ha('right')
                text.set_x(shifts[i])
            ax.set_aspect(1.0 / ax.get_data_ratio())
            fig.tight_layout(pad=0.5, w_pad=0, h_pad=0)
            for fmt in args.file_format:
                fig.savefig('{}/{}_roc_auc.{}'.format(args.out_dir, model_name,
                                                      fmt),
                            format=fmt, bbox_inches='tight',
                            # matplotlib GH#15497
                            dpi='figure' if fmt == 'pdf' else fig_dpi)
            pd.DataFrame(tsv_scores).to_csv(
                '{}/{}_roc_auc.tsv'.format(args.out_dir, model_name),
                index=False, sep='\t')

            # pr curves
            if data_type == 'kraken':
                colors = ['purplish', 'steel grey']
            elif data_type == 'htseq':
                colors = ['turquoise', 'steel grey']
            else:
                colors = ['purplish', 'burnt orange', 'dark sky blue']

            colors = sns.xkcd_palette(colors)

            tsv_scores = {k: [] for k in ['dtype', 'split', 'rec', 'pre']}
            fig, ax = plt.subplots(figsize=(fig_dim, fig_dim))
            for ridx, _ in enumerate(split_results):
                pres, pr_scores = [], []
                mean_recs = np.linspace(0, 1, 1000)
                for split_idx, split_result in enumerate(split_results[ridx]):
                    if split_result is None:
                        continue
                    rec = split_result['scores']['te']['rec'][::-1]
                    pre = split_result['scores']['te']['pre'][::-1]
                    pres.append(np.interp(mean_recs, rec, pre))
                    pr_scores.append(split_result['scores']['te']['pr_auc'])
                    if ridx == 0:
                        tsv_data_type = data_type
                    elif data_type == 'combo':
                        tsv_data_type = 'htseq' if ridx == 1 else 'kraken'
                    else:
                        tsv_data_type = 'clinical'
                    tsv_scores['dtype'].extend([tsv_data_type] * len(rec))
                    tsv_scores['split'].extend([split_idx + 1] * len(rec))
                    tsv_scores['rec'].extend(rec)
                    tsv_scores['pre'].extend(pre)
                mean_pres = np.mean(pres, axis=0)
                std_pres = np.std(pres, axis=0)
                pres_upper = np.minimum(mean_pres + std_pres, 1)
                pres_lower = np.maximum(mean_pres - std_pres, 0)
                if data_type == 'combo':
                    label = '+'.join([dtype_labels[ridx], dtype_labels[-1]])
                    zorder = 2.5 if ridx == 0 else 2.2 if ridx == 1 else 2
                elif ridx == 0:
                    label = '+'.join([dtype_labels[ridx], dtype_labels[-1]])
                    zorder = 2.5
                else:
                    label = dtype_labels[-1]
                    zorder = 2
                ax.step(mean_recs, mean_pres, alpha=0.8, color=colors[ridx],
                        label=('AUPRC = {:.2f}'.format(np.mean(pr_scores))
                               if ridx == 0 else None), lw=2, where='post',
                        zorder=zorder)
                ax.fill_between(mean_recs, pres_lower, pres_upper, alpha=0.1,
                                color=colors[ridx], zorder=zorder)
            ax.set_title(figure_title, loc='center', pad=8,
                         fontdict={'fontsize': title_fontsize})
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
            ax.tick_params(which='major', length=5, width=1)
            ax.tick_params(which='minor', width=1)
            ax.margins(0)
            ax.grid(False)
            legend = ax.legend(loc='lower right', borderpad=0.1,
                               borderaxespad=0.1, frameon=False,
                               labelspacing=0.2, fontsize=legend_fontsize)
            legend.set_title(
                '+'.join([abbr_dtype_labels[0], abbr_dtype_labels[-1]]),
                prop={'weight': 'regular', 'size': legend_fontsize})
            legend._legend_box.align = 'right'
            for item in legend.legendHandles:
                item.set_visible(False)
            text_widths = [
                text.get_window_extent(fig.canvas.get_renderer()).width
                for text in legend.get_texts()]
            max_width = max(text_widths)
            shifts = [max_width - w for w in text_widths]
            for i, text in enumerate(legend.get_texts()):
                text.set_ha('right')
                text.set_x(shifts[i])
            ax.set_aspect(1.0 / ax.get_data_ratio())
            fig.tight_layout(pad=0.5, w_pad=0, h_pad=0)
            for fmt in args.file_format:
                fig.savefig('{}/{}_pr_auc.{}'.format(args.out_dir, model_name,
                                                     fmt),
                            format=fmt, bbox_inches='tight',
                            # matplotlib GH#15497
                            dpi='figure' if fmt == 'pdf' else fig_dpi)
            pd.DataFrame(tsv_scores).to_csv(
                '{}/{}_pr_auc.tsv'.format(args.out_dir, model_name),
                index=False, sep='\t')

            # num selected features vs scores
            param_cv_scores = load('{}/resp/{name}/{name}_param_cv_scores.pkl'
                                   .format(model_results_dir, name=model_name))

            if data_type == 'kraken':
                colors = ['dark sky blue', 'purplish']
            elif data_type == 'htseq':
                colors = ['burnt orange', 'turquoise']
            else:
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
                fig, ax = plt.subplots(figsize=(fig_dim, fig_dim))
                if param_parts[-1] in ('k', 'n_features_to_select'):
                    param_ext = 'k'
                    x_label = 'Num selected features'
                    x_axis = np.insert(np.linspace(2, 400, num=200, dtype=int),
                                       0, 1)
                    ax.get_xaxis().set_major_locator(ticker.FixedLocator(
                        [1, 100, 200, 300, 400]))
                    ax.get_xaxis().set_minor_locator(ticker.FixedLocator(
                        [50, 150, 250, 350]))
                    ax.grid(True, alpha=0.3, which='both')
                elif param_parts[-1] == 'C':
                    param_ext = 'c'
                    x_label = 'C'
                    x_axis = (np.logspace(-2, 3, 6) if data_type == 'kraken'
                              else np.logspace(-2, 1, 4))
                    ax.set_xscale('log')
                    ax.set_xticks(x_axis)
                    ax.get_xaxis().set_minor_locator(ticker.LogLocator(
                        base=10, subs='all', numticks=8))
                    ax.get_xaxis().set_major_locator(ticker.LogLocator(
                        base=10, numticks=len(x_axis)))
                    ax.grid(True, alpha=0.3, which='major')
                elif param_parts[-1] == 'l1_ratio':
                    param_ext = 'l1r'
                    x_label = 'L1 ratio'
                    x_axis = np.array([0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95,
                                       0.99, 1.])
                    ax.set_xticks([0.1, 0.3, 0.5, 0.7, 0.9, 1])
                    ax.get_xaxis().set_major_formatter(ticker.FixedFormatter(
                        ['0.1', '0.3', '0.5', '0.7', '0.9', '1']))
                    ax.grid(True, alpha=0.3, which='major')
                l_metric_labels = [s.lower() for s in metric_labels]
                tsv_scores = {k: [] for k in [param_ext] + l_metric_labels}
                for metric_idx, metric in enumerate(metrics):
                    mean_cv_scores, std_cv_scores = [], []
                    param_metric_scores = (
                        param_cv_scores[param][metric]['scores'])
                    param_metric_stdev = (
                        param_cv_scores[param][metric]['stdev'])
                    for idx, param_value_scores in enumerate(
                            param_metric_scores):
                        if metric_idx == 0:
                            tsv_scores[param_ext].extend(
                                [x_axis[idx]] * len(param_value_scores))
                        tsv_scores[l_metric_labels[metric_idx]].extend(
                            param_value_scores)
                        mean_cv_scores.append(np.mean(param_value_scores))
                        std_cv_scores.append(np.std(param_value_scores))
                    zorder = (2.5 if metric_idx == 0 else
                              2.2 if metric_idx == 1 else 2)
                    ax.plot(x_axis, mean_cv_scores,
                            color=colors[metric_idx], lw=2, alpha=0.8,
                            label='{}'.format(metric_labels[metric_idx]),
                            zorder=zorder)
                    ax.fill_between(
                        x_axis,
                        [m - s for m, s in zip(mean_cv_scores, std_cv_scores)],
                        [m + s for m, s in zip(mean_cv_scores, std_cv_scores)],
                        alpha=0.1, color=colors[metric_idx], zorder=zorder)
                ax.set_title(figure_title, loc='center', pad=8,
                             fontdict={'fontsize': title_fontsize})
                ax.set_xlabel(x_label, fontsize=axis_fontsize)
                ax.set_ylabel('Score', fontsize=axis_fontsize)
                ax.set_xlim([min(x_axis), max(x_axis)])
                ax.set_ylim([0.0, 1.0])
                ax.set_yticks(np.arange(0.0, 1.1, 0.2))
                ax.get_yaxis().set_major_formatter(ticker.FixedFormatter(
                    ['0', '0.2', '0.4', '0.6', '0.8', '1']))
                ax.tick_params(axis='both', labelsize=axis_fontsize)
                ax.tick_params(which='major', length=5, width=1)
                ax.tick_params(which='minor', length=3, width=1)
                ax.margins(0)
                legend = ax.legend(loc='lower right', borderpad=0.2,
                                   borderaxespad=0.1, fontsize=legend_fontsize,
                                   handlelength=1.5, handletextpad=0.3,
                                   labelspacing=0.2)
                legend._legend_box.align = 'right'
                # text_widths = [
                #     text.get_window_extent(fig.canvas.get_renderer()).width
                #     for text in legend.get_texts()]
                # max_width = max(text_widths)
                # shifts = [max_width - w for w in text_widths]
                # for i, text in enumerate(legend.get_texts()):
                #     text.set_ha('right')
                #     text.set_x(shifts[i])
                ax.set_aspect(1.0 / ax.get_data_ratio())
                fig.tight_layout(pad=0.5, w_pad=0, h_pad=0)
                for fmt in args.file_format:
                    fig.savefig('{}/{}_{}_vs_score.{}'.format(
                        args.out_dir, model_name, param_ext, fmt),
                                format=fmt, bbox_inches='tight',
                                # matplotlib GH#15497
                                dpi='figure' if fmt == 'pdf' else fig_dpi)
                pd.DataFrame(tsv_scores).to_csv(
                    '{}/{}_{}_vs_score.tsv'.format(args.out_dir, model_name,
                                                   param_ext),
                    index=False, sep='\t')
