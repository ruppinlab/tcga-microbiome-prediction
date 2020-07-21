import os
import warnings
from argparse import ArgumentParser
from glob import glob

warnings.filterwarnings('ignore', category=FutureWarning,
                        module='sklearn.utils.deprecation')
warnings.filterwarnings('ignore', category=FutureWarning,
                        module='rpy2.robjects.pandas2ri')

import numpy as np
import pandas as pd
import rpy2.rinterface_lib.embedded as r_embedded

r_embedded.set_initoptions(
    ('rpy2', '--quiet', '--no-save', '--max-ppsize=500000'))

import rpy2.robjects as robjects
from joblib import delayed, dump, Parallel
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder, StandardScaler
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.util import Surv
from tabulate import tabulate

from sksurv_extensions.model_selection import (
    SurvivalStratifiedShuffleSplit,
    SurvivalStratifiedSampleFromGroupShuffleSplit)

numpy2ri.activate()
pandas2ri.activate()


def fit_models(X, y, groups, group_weights):
    test_splits = 100
    test_size = 0.25
    random_seed = 777

    pipe = Pipeline([('trf0', StandardScaler()),
                     ('srv1', CoxPHSurvivalAnalysis(alpha=1e-09, n_iter=10000,
                                                    ties='efron', tol=1e-09))])

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

    split_models, split_results = [], []
    for train_idxs, test_idxs in cv.split(X, y, groups, **cv_split_params):
        try:
            pipe.fit(X.iloc[train_idxs], y[train_idxs])
            score = pipe.score(X.iloc[test_idxs], y[test_idxs])
            y_pred = pipe.predict(X.iloc[test_idxs])
            surv_funcs = pipe.predict_survival_function(X.iloc[test_idxs])
        except Exception:
            split_models.append(None)
            split_results.append(None)
        else:
            split_models.append(pipe)
            split_scores = {'te': {'score': score, 'y_pred': y_pred}}
            split_results.append({'scores': split_scores,
                                  'surv_funcs': surv_funcs})

    return split_models, split_results


r_base = importr('base')
r_biobase = importr('Biobase')

ordinal_encoder_categories = {
    'tumor_stage': ['NA', 'x', 'i', 'i or ii', 'ii', 'iii', 'iv']}

parser = ArgumentParser()
parser.add_argument('--data-dir', type=str, default='data', help='data dir')
parser.add_argument('--out-dir', type=str, default=os.getcwd(), help='out dir')
parser.add_argument('--n-jobs', type=int, default=-1, help='num parallel jobs')
parser.add_argument('--verbose', type=int, default=1, help='verbosity')
args = parser.parse_args()

all_X, all_y, all_groups, all_group_weights = [], [], [], []
eset_files = sorted(glob('{}/tcga_*_surv_*_eset.rds'.format(args.data_dir)))
num_esets = len(eset_files)
for eset_idx, eset_file in enumerate(eset_files):
    file_basename = os.path.splitext(os.path.split(eset_file)[1])[0]

    if args.verbose < 2:
        print('Loading {:d}/{:d} esets'.format(eset_idx + 1, num_esets),
              end='\r', flush=True)
    else:
        print('Loading {}'.format(file_basename))

    eset = r_base.readRDS(eset_file)
    sample_meta = r_biobase.pData(eset)
    X = pd.DataFrame(index=sample_meta.index)
    y = Surv.from_dataframe('Status', 'Survival_in_days', sample_meta)

    if 'Group' in sample_meta.columns:
        groups = np.array(sample_meta['Group'], dtype=int)
        if ('GroupWeight' in sample_meta.columns
                and sample_meta['GroupWeight'].unique().size > 1):
            group_weights = np.array(sample_meta['GroupWeight'], dtype=float)
        else:
            group_weights = None
    else:
        groups = None
        group_weights = None

    if sample_meta['gender'].unique().size > 1:
        ohe = OneHotEncoder(drop='first', sparse=False)
        ohe.fit(sample_meta[['gender']])
        feature_name = 'gender_{}'.format(ohe.categories_[0][1])
        X[feature_name] = ohe.transform(sample_meta[['gender']])
    X['age_at_diagnosis'] = sample_meta[['age_at_diagnosis']]
    if sample_meta['tumor_stage'].unique().size > 1:
        ode = OrdinalEncoder(categories=[
            ordinal_encoder_categories['tumor_stage']])
        ode.fit(sample_meta[['tumor_stage']])
        X['tumor_stage'] = ode.transform(sample_meta[['tumor_stage']])

    all_X.append(X)
    all_y.append(y)
    all_groups.append(groups)
    all_group_weights.append(group_weights)

if args.verbose < 2:
    print(flush=True)

print('Running Cox models')
all_models, all_results = zip(*Parallel(
    n_jobs=args.n_jobs, verbose=args.verbose)(
        delayed(fit_models)(X, y, groups, group_weights)
        for X, y, groups, group_weights in
        zip(all_X, all_y, all_groups, all_group_weights)))

if args.verbose < 1:
    print(flush=True)

mean_scores = []
all_scores_df = None
for eset_file, split_models, split_results in zip(eset_files, all_models,
                                                  all_results):
    file_basename = os.path.splitext(os.path.split(eset_file)[1])[0]
    _, cancer, analysis, target, data_type, *rest = file_basename.split('_')
    data_type = 'expr' if data_type == 'htseq' else data_type

    scores = []
    for split_result in split_results:
        if split_result is not None:
            scores.append(split_result['scores']['te']['score'])
        else:
            scores.append(np.nan)

    mean_score = np.nanmean(scores)
    mean_scores.append([analysis, cancer, target, data_type, mean_score])

    dataset_name = '_'.join(file_basename.split('_')[:-1])
    model_name = '_'.join([dataset_name, 'cox'])
    dump(split_models, '{}/{}_split_models.pkl'.format(args.out_dir,
                                                       model_name))
    dump(split_results, '{}/{}_split_results.pkl'.format(args.out_dir,
                                                         model_name))

    scores_df = pd.DataFrame({dataset_name: scores})
    if all_scores_df is None:
        all_scores_df = scores_df
    else:
        all_scores_df = pd.concat([all_scores_df, scores_df], axis=1)

all_scores_df.to_csv('{}/cox_covariate_model_scores.tsv'.format(args.out_dir),
                     sep='\t')

dump(all_scores_df, '{}/cox_covariate_model_scores.pkl'.format(args.out_dir))

r_base.saveRDS(all_scores_df,
               '{}/cox_covariate_model_scores.rds'.format(args.out_dir))

mean_scores_df = pd.DataFrame(mean_scores, columns=[
    'Analysis', 'Cancer', 'Target', 'Data Type', 'Mean Score'])
mean_scores_df.to_csv('{}/cox_covariate_model_mean_scores.tsv'
                      .format(args.out_dir), index=False, sep='\t')
if args.verbose > 0:
    print(tabulate(
        mean_scores_df.sort_values(['Analysis', 'Cancer', 'Target',
                                    'Data Type']),
        floatfmt='.4f', showindex=False, headers='keys'))
