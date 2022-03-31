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
import rpy2.rinterface_lib.embedded as r_embedded

r_embedded.set_initoptions(
    ('rpy2', '--quiet', '--no-save', '--max-ppsize=500000'))

import rpy2.robjects as robjects
import seaborn as sns
from joblib import delayed, load, Parallel
from matplotlib import ticker
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr

from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder
from sksurv.metrics import cumulative_dynamic_auc
from sksurv.util import Surv

from sksurv_extensions.model_selection import (
    SurvivalStratifiedShuffleSplit,
    SurvivalStratifiedSampleFromGroupShuffleSplit)

numpy2ri.activate()
pandas2ri.activate()


def get_eset_dataset(eset_file):
    eset = r_base.readRDS(eset_file)
    sample_meta = r_biobase.pData(eset)
    X = pd.DataFrame(index=sample_meta.index)
    y = Surv.from_dataframe(sample_meta_stat_col, sample_meta_surv_col,
                            sample_meta)

    if 'Group' in sample_meta.columns:
        groups = np.array(sample_meta['Group'], dtype=int)
        if ('GroupWeight' in sample_meta.columns
                and sample_meta['GroupWeight'].unique().size > 1):
            group_weights = np.array(sample_meta['GroupWeight'],
                                     dtype=float)
        else:
            group_weights = None
    else:
        groups = None
        group_weights = None

    X['age_at_diagnosis'] = sample_meta[['age_at_diagnosis']]
    if sample_meta['gender'].unique().size > 1:
        ohe = OneHotEncoder(drop='first', sparse=False)
        ohe.fit(sample_meta[['gender']])
        feature_name = 'gender_{}'.format(ohe.categories_[0][1])
        X[feature_name] = ohe.transform(sample_meta[['gender']])
    if sample_meta['tumor_stage'].unique().size > 1:
        ode = OrdinalEncoder(categories=[
            ordinal_encoder_categories['tumor_stage']])
        ode.fit(sample_meta[['tumor_stage']])
        X['tumor_stage'] = ode.transform(
            sample_meta[['tumor_stage']])

    return X, y, groups, group_weights


def get_cv_split_idxs(X, y, groups, group_weights):
    if groups is None:
        cv = SurvivalStratifiedShuffleSplit(
            n_splits=test_splits, test_size=test_size,
            random_state=random_seed)
        cv_split_params = {}
    else:
        cv = SurvivalStratifiedSampleFromGroupShuffleSplit(
            n_splits=test_splits, test_size=test_size,
            random_state=random_seed)
        cv_split_params = {'weights': group_weights}

    cv_split_idxs = []
    for train_idxs, test_idxs in cv.split(X, y, groups, **cv_split_params):
        cv_split_idxs.append((train_idxs, test_idxs))

    return cv_split_idxs


warnings.filterwarnings('ignore', category=RuntimeWarning,
                        message='^overflow encountered in power',
                        module='sksurv.linear_model.coxph')
warnings.filterwarnings('ignore', category=RuntimeWarning,
                        message='^invalid value encountered in true_divide',
                        module='sksurv.metrics')

# suppress linux conda qt5 wayland warning
if sys.platform.startswith('linux'):
    os.environ['XDG_SESSION_TYPE'] = 'x11'

parser = ArgumentParser()
parser.add_argument('--data-dir', type=str, default='data')
parser.add_argument('--results-dir', type=str, default='results',
                    help='results dir')
parser.add_argument('--out-dir', type=str, default='figures/td_auc',
                    help='out dir')
parser.add_argument('--test-splits', type=int, help='num test splits')
parser.add_argument('--test-size', type=float, help='test split size')
parser.add_argument('--file-format', type=str, nargs='+',
                    choices=['png', 'pdf', 'svg', 'tif'], default=['png'],
                    help='save file format')
parser.add_argument('--n-jobs', type=int, default=-1, help='num parallel jobs')
parser.add_argument('--verbose', type=int, default=0, help='verbosity')
args = parser.parse_args()

model_results_dir = '{}/models'.format(args.results_dir)

os.makedirs(args.out_dir, mode=0o755, exist_ok=True)

test_splits = 100 if args.test_splits is None else args.test_splits
test_size = 0.25 if args.test_size is None else args.test_size
random_seed = 777

ordinal_encoder_categories = {
    'tumor_stage': ['0', 'i', 'i or ii', 'ii', 'NA', 'iii', 'iv']}
sample_meta_stat_col = 'Status'
sample_meta_surv_col = 'Survival_in_days'

title_fontsize = 14
axis_fontsize = 12
legend_fontsize = 12
fig_let_fontsize = 48
fig_dim = 4
fig_dpi = 300
time_interval_days = 30
days_per_year = 365.2422

plt.rcParams['figure.max_open_warning'] = 0
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'sans']

r_base = importr('base')
r_biobase = importr('Biobase')

split_results_regex = re.compile('^(.+?_cnet)_split_results\\.pkl$')
for dirpath, dirnames, filenames in sorted(os.walk(model_results_dir)):
    for filename in filenames:
        if m := re.search(split_results_regex, filename):
            model_name = m.group(1)
            print(model_name)
            _, cancer, analysis, target, data_type, *rest = (
                model_name.split('_'))
            figure_title = '{} {}'.format(cancer.upper(), target.upper())
            data_type_label = ('Expression' if data_type == 'htseq' else
                               'Microbiome')

            eset_files, split_results = [], []
            dataset_name = '_'.join(model_name.split('_')[:-1])
            eset_files.append('{}/{}_eset.rds'.format(
                args.data_dir, dataset_name))
            split_results.append(load(
                '{}/surv/{name}/{name}_split_results.pkl'
                .format(model_results_dir, name=model_name)))
            if data_type in ('kraken', 'htseq'):
                cox_model_name = '_'.join([dataset_name, 'cox', 'clinical'])
                eset_files.append('{}/{}_eset.rds'
                                  .format(args.data_dir, dataset_name))
                split_results.append(load(
                    '{}/surv/{name}/{name}_split_results.pkl'
                    .format(model_results_dir, name=cox_model_name)))
            else:
                for new_data_type in ('htseq_counts', 'kraken'):
                    new_model_name = '_'.join(
                        model_name.split('_')[:-2] + [new_data_type, 'cnet'])
                    new_dataset_name = '_'.join(new_model_name.split('_')[:-1])
                    eset_files.append('{}/{}_eset.rds'.format(
                        args.data_dir, new_dataset_name))
                    split_results.append(load(
                        '{}/surv/{name}/{name}_split_results.pkl'
                        .format(model_results_dir, name=new_model_name)))

            datasets = [get_eset_dataset(file) for file in eset_files]
            all_cv_split_idxs = Parallel(
                n_jobs=args.n_jobs, verbose=args.verbose)(
                    delayed(get_cv_split_idxs)(X, y, groups, group_weights)
                    for X, y, groups, group_weights in datasets)

            if data_type == 'kraken':
                colors = ['dark sky blue', 'purplish']
            elif data_type == 'htseq':
                colors = ['burnt orange']
            else:
                colors = ['purplish', 'burnt orange', 'dark sky blue']
            colors.append('steel grey')
            colors = sns.xkcd_palette(colors)

            # time-dependent cumulative/dynamic AUCs
            tsv_scores = {k: [] for k in ['data_type', 'time', 'auc']}
            fig, ax = plt.subplots(figsize=(fig_dim, fig_dim), dpi=fig_dpi)
            for ridx, _ in enumerate(split_results):
                y = datasets[ridx][1]
                y_stat, y_time = y.dtype.names
                times, aucs = [], []
                for split_idx, (train_idxs, test_idxs) in enumerate(
                        all_cv_split_idxs[ridx]):
                    split_result = split_results[ridx][split_idx]
                    if split_result is None:
                        continue
                    time_test_idxs = test_idxs[
                        (y[test_idxs][y_time]
                         >= np.min(y[train_idxs][y_time]))
                        & (y[test_idxs][y_time]
                           <= np.max(y[train_idxs][y_time]))]
                    if not any(y[time_test_idxs][y_stat]):
                        print('Split {} has all censored y_test samples '
                              'within y_train times, skipping'
                              .format(split_idx))
                        continue
                    y_pred_idxs = np.where(np.isin(test_idxs, time_test_idxs,
                                                   assume_unique=True))[0]
                    y_pred = (
                        split_result['scores']['te']['y_pred'][y_pred_idxs])
                    time = np.arange(np.min(y[time_test_idxs][y_time]),
                                     np.max(y[time_test_idxs][y_time]),
                                     time_interval_days)
                    auc = cumulative_dynamic_auc(
                        y[train_idxs], y[time_test_idxs], y_pred, time)[0]
                    nan_auc = np.isnan(auc)
                    if np.all(nan_auc):
                        continue
                    time = time[np.logical_not(nan_auc)]
                    auc = auc[np.logical_not(nan_auc)]
                    times.append(time)
                    aucs.append(auc)
                    if ridx == 0:
                        tsv_data_type = data_type
                    elif data_type == 'combo':
                        tsv_data_type = 'htseq' if ridx == 1 else 'kraken'
                    else:
                        tsv_data_type = 'clinical'
                    tsv_scores['data_type'].extend([tsv_data_type] * len(time))
                    tsv_scores['time'].extend(time)
                    tsv_scores['auc'].extend(auc)
                interp_aucs = []
                mean_times = np.linspace(
                    min(t[0] for t in times), max(t[-1] for t in times), 1000)
                for time, auc in zip(times, aucs):
                    interp_aucs.append(np.interp(mean_times, time, auc))
                mean_times = mean_times / days_per_year
                mean_aucs = np.mean(interp_aucs, axis=0)
                std_aucs = np.std(interp_aucs, axis=0)
                aucs_upper = np.minimum(mean_aucs + std_aucs, 1)
                aucs_lower = np.maximum(mean_aucs - std_aucs, 0)
                if data_type == 'combo':
                    dtype_label = ('Combo' if ridx == 0 else
                                   'Expression' if ridx == 1 else 'Microbiome')
                    label = '{} + Clinical'.format(dtype_label)
                    color = colors[ridx]
                    zorder = 2.5 if ridx == 0 else 2.2 if ridx == 1 else 2
                elif ridx == 0:
                    label = '{} + Clinical'.format(data_type_label)
                    color = colors[0]
                    zorder = 2.5
                else:
                    label = 'Clinical'
                    color = colors[-1]
                    zorder = 2
                ax.plot(mean_times, mean_aucs, alpha=0.8, color=color, lw=2,
                        label=r'{} AUC = $\bf{{{:.2f}}}$'.format(
                            label, np.mean(mean_aucs)), zorder=zorder)
                ax.fill_between(mean_times, aucs_lower, aucs_upper, alpha=0.1,
                                color=color, zorder=zorder)
                ax.axhline(np.mean(mean_aucs), alpha=0.5, color=color,
                           linestyle='--', lw=1.5, zorder=1)
            xaxis_tick_base = (3 if max(mean_times) > 20 else
                               2 if max(mean_times) > 10 else 1)
            ax.set_title(figure_title, loc='left', y=1.0, pad=4,
                         fontdict={'fontsize': title_fontsize})
            ax.set_xlabel('Years from diagnosis', fontsize=axis_fontsize,
                          labelpad=5)
            ax.get_xaxis().set_major_locator(
                ticker.MultipleLocator(base=xaxis_tick_base))
            ax.set_ylabel('Time-dependent AUC', fontsize=axis_fontsize,
                          labelpad=5)
            ax.set_yticks(np.arange(0.2, 1.1, 0.2))
            ax.get_yaxis().set_major_formatter(ticker.FixedFormatter(
                ['0.2', '0.4', '0.6', '0.8', '1']))
            ax.set_ylim([0.2, 1.01])
            ax.tick_params(axis='both', labelsize=axis_fontsize)
            ax.tick_params(which='major', width=1)
            ax.tick_params(which='major', length=5)
            ax.tick_params(which='minor', width=1)
            ax.margins(0)
            ax.grid(False)
            legend = ax.legend(loc='lower right', frameon=False, borderpad=0.1,
                               prop={'size': legend_fontsize})
            # legend.set_title(figure_title, prop={'weight': 'regular',
            #                                      'size': axis_fontsize})
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
                fig.savefig('{}/{}_td_auc.{}'.format(args.out_dir, model_name,
                                                     fmt),
                            format=fmt, bbox_inches='tight')
            pd.DataFrame(tsv_scores).to_csv(
                '{}/{}_td_auc.tsv'.format(args.out_dir, model_name),
                index=False, sep='\t')
