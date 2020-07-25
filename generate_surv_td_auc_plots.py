import os
import re
import sys
import warnings
from argparse import ArgumentParser

warnings.filterwarnings('ignore', category=FutureWarning,
                        module='sklearn.utils.deprecation')
warnings.filterwarnings('ignore', category=FutureWarning,
                        module='rpy2.robjects.pandas2ri')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.api.types import (is_categorical_dtype, is_object_dtype,
                              is_string_dtype)
import rpy2.rinterface_lib.embedded as r_embedded

r_embedded.set_initoptions(
    ('rpy2', '--quiet', '--no-save', '--max-ppsize=500000'))

import rpy2.robjects as robjects
import seaborn as sns
from joblib import load, Memory
from matplotlib import ticker
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr

from sklearn.base import BaseEstimator, clone
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder
from sksurv.metrics import (cumulative_dynamic_auc, brier_score,
                            integrated_brier_score)
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.util import Surv

from sksurv_extensions.model_selection import (
    SurvivalStratifiedShuffleSplit,
    SurvivalStratifiedSampleFromGroupShuffleSplit)

numpy2ri.activate()
pandas2ri.activate()

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
parser.add_argument('--out-dir', type=str, default='figures', help='out dir')
parser.add_argument('--test-splits', type=int, help='num test splits')
parser.add_argument('--test-size', type=float, help='test split size')
parser.add_argument('--file-format', type=str, nargs='+',
                    choices=['png', 'pdf', 'svg', 'tif'], default=['png'],
                    help='save file format')
args = parser.parse_args()

os.makedirs(args.out_dir, mode=0o755, exist_ok=True)

test_splits = 100 if args.test_splits is None else args.test_splits
test_size = 0.25 if args.test_size is None else args.test_size
random_seed = 777

r_base = importr('base')
r_biobase = importr('Biobase')

penalty_factor_meta_col = 'Penalty Factor'
sample_meta_cols = ['age_at_diagnosis', 'gender', 'tumor_stage']
ordinal_encode_cols = ['tumor_stage']
ordinal_encoder_categories = {
    'tumor_stage': ['NA', 'x', 'i', 'i or ii', 'ii', 'iii', 'iv']}

title_fontsize = 16
axis_fontsize = 12
legend_fontsize = 10
fig_let_fontsize = 48
fig_dim = 4
fig_dpi = 300
time_interval_days = 30
days_per_year = 365.2422

plt.rcParams['figure.max_open_warning'] = 0
plt.rcParams['font.family'] = ['Nimbus Sans']

fig_count = {}
results_dirname_regex = re.compile('^(tcga_.+?_cnet)$')
for dirpath, dirnames, filenames in sorted(os.walk(args.results_dir)):
    for dirname in dirnames:
        if m := re.search(results_dirname_regex, dirname):
            model_name = m.group(1)
            print(model_name)
            _, cancer, analysis, target, data_type, *rest = (
                model_name.split('_'))
            legend_title = '{} {}'.format(cancer.upper(), target.upper())
            data_type_label = ('Expression' if data_type == 'htseq' else
                               'Microbiome')

            split_results = []
            split_results.append(load(
                '{}/surv/{name}/{name}_split_results.pkl'
                .format(args.results_dir, name=model_name)))
            dataset_name = '_'.join(model_name.split('_')[:-1])
            cox_model_name = '_'.join([dataset_name, 'cox'])
            split_results.append(load(
                '{}/surv/{name}/{name}_split_results.pkl'
                .format(args.results_dir, name=cox_model_name)))

            if data_type == 'kraken':
                fig_num = '1'
                colors = ['dark sky blue']
            else:
                # fig_num Ex1 or Ex2 below depending on OS or PFI
                colors = ['burnt orange']
            colors.append('steel grey')
            colors = sns.xkcd_palette(colors)

            eset_file = '{}/{}_eset.rds'.format(args.data_dir, dataset_name)
            eset = r_base.readRDS(eset_file)
            sample_meta = r_biobase.pData(eset)
            X = pd.DataFrame(index=sample_meta.index)
            y = Surv.from_dataframe('Status', 'Survival_in_days', sample_meta)

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
                X['tumor_stage'] = ode.transform(sample_meta[['tumor_stage']])

            y_stat = y.dtype.names[0]
            y_time = y.dtype.names[1]

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
            for train_idxs, test_idxs in cv.split(X, y, groups,
                                                  **cv_split_params):
                cv_split_idxs.append((train_idxs, test_idxs))

            # time-dependent AUCs
            fig, ax = plt.subplots(figsize=(fig_dim, fig_dim), dpi=fig_dpi)
            for ridx, _ in enumerate(split_results):
                times, aucs = [], []
                for split_idx, (train_idxs, test_idxs) in enumerate(
                        cv_split_idxs):
                    split_result = split_results[ridx][split_idx]
                    if split_result is None:
                        continue
                    time_test_idxs = test_idxs[
                        (y[test_idxs][y_time]
                         >= np.min(y[train_idxs][y_time]))
                        & (y[test_idxs][y_time]
                           <= np.max(y[train_idxs][y_time]))]
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
                    times.append(time[np.logical_not(nan_auc)])
                    aucs.append(auc[np.logical_not(nan_auc)])
                interp_aucs = []
                mean_times = np.linspace(
                    min(t[0] for t in times), max(t[-1] for t in times), 1000)
                for time, auc in zip(times, aucs):
                    interp_aucs.append(np.interp(mean_times, time, auc))
                mean_times = mean_times / days_per_year
                mean_auc = np.mean(interp_aucs, axis=0)
                std_auc = np.std(interp_aucs, axis=0)
                aucs_upper = np.minimum(mean_auc + std_auc, 1)
                aucs_lower = np.maximum(mean_auc - std_auc, 0)
                if ridx == 0:
                    label = '{} + Covariate'.format(data_type_label)
                    color = colors[0]
                    zorder = 2.5
                else:
                    label = 'Covariate'
                    color = colors[-1]
                    zorder = 2
                ax.plot(mean_times, mean_auc, alpha=0.8, color=color, lw=2,
                        label=r'{} AUC = $\bf{{{:.2f}}}$'.format(
                            label, np.mean(mean_auc)), zorder=zorder)
                ax.fill_between(mean_times, aucs_lower, aucs_upper, alpha=0.1,
                                color=color, zorder=zorder)
                ax.axhline(np.mean(mean_auc), alpha=0.5, color=color,
                           linestyle='--', lw=1.5, zorder=1)
            xaxis_tick_base = (3 if max(mean_times) > 20 else
                               2 if max(mean_times) > 10 else 1)
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
            legend.set_title(legend_title, prop={'weight': 'bold',
                                                 'size': axis_fontsize})
            legend._legend_box.align = 'right'
            for item in legend.legendHandles:
                item.set_visible(False)
            renderer = fig.canvas.get_renderer()
            shift = max([text.get_window_extent(renderer).width
                         for text in legend.get_texts()])
            for text in legend.get_texts():
                text.set_ha('right')
                text.set_position((shift, 0))
            ax.set_aspect(1.0 / ax.get_data_ratio())
            fig.tight_layout(pad=0.5, w_pad=0, h_pad=0)
            if data_type == 'htseq':
                fig_num = 'Ex1' if target == 'os' else 'Ex2'
            fig_num = '{}C'.format(fig_num)
            if fig_num not in fig_count:
                fig_count[fig_num] = 1
            for fmt in args.file_format:
                fig.savefig('{}/Figure_{}{:02d}.{}'.format(
                    args.out_dir, fig_num, fig_count[fig_num], fmt),
                            format=fmt, bbox_inches='tight')
            fig_count[fig_num] += 1
