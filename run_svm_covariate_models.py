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

from sklearn.metrics import (
    auc, average_precision_score, balanced_accuracy_score,
    precision_recall_curve, roc_auc_score, roc_curve)
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder, StandardScaler
from sklearn.svm import SVC
from tabulate import tabulate

from sklearn_extensions.model_selection import RepeatedStratifiedGroupKFold

numpy2ri.activate()
pandas2ri.activate()


def calculate_test_scores(pipe, X_test, y_test, pipe_predict_params,
                          test_sample_weights=None):
    scores = {}
    if hasattr(pipe, 'decision_function'):
        y_score = pipe.decision_function(X_test, **pipe_predict_params)
    else:
        y_score = pipe.predict_proba(X_test, **pipe_predict_params)[:, 1]
    scores['y_score'] = y_score
    for metric in metrics:
        if metric == 'roc_auc':
            scores[metric] = roc_auc_score(
                y_test, y_score, sample_weight=test_sample_weights)
            scores['fpr'], scores['tpr'], _ = roc_curve(
                y_test, y_score, pos_label=1,
                sample_weight=test_sample_weights)
        elif metric == 'balanced_accuracy':
            y_pred = pipe.predict(X_test, **pipe_predict_params)
            scores['y_pred'] = y_pred
            scores[metric] = balanced_accuracy_score(
                y_test, y_pred, sample_weight=test_sample_weights)
        elif metric == 'average_precision':
            scores[metric] = average_precision_score(
                y_test, y_score, sample_weight=test_sample_weights)
            scores['pre'], scores['rec'], _ = precision_recall_curve(
                y_test, y_score, pos_label=1,
                sample_weight=test_sample_weights)
            scores['pr_auc'] = auc(scores['rec'], scores['pre'])
    return scores


def fit_models(X, y, groups, sample_weights, cv_split_flag):
    random_seed = 777

    pipe = Pipeline([('trf0', StandardScaler()),
                     ('clf1', SVC(kernel='linear', class_weight='balanced',
                                  random_state=random_seed))])

    if cv_split_flag:
        test_splits = 3
        test_repeats = 33
    else:
        test_splits = 4
        test_repeats = 25

    if groups is None:
        cv = RepeatedStratifiedKFold(n_splits=test_splits,
                                     n_repeats=test_repeats,
                                     random_state=random_seed)
    else:
        cv = RepeatedStratifiedGroupKFold(n_splits=test_splits,
                                          n_repeats=test_repeats,
                                          random_state=random_seed)

    split_results = []
    for train_idxs, test_idxs in cv.split(X, y, groups):
        train_sample_weights = None
        test_sample_weights = None
        if sample_weights is not None:
            train_sample_weights = sample_weights[train_idxs]
            test_sample_weights = sample_weights[test_idxs]
        pipe.fit(X.iloc[train_idxs], y[train_idxs],
                 clf1__sample_weight=train_sample_weights)
        split_scores = {'te': calculate_test_scores(
            pipe, X.iloc[test_idxs], y[test_idxs], pipe_predict_params={},
            test_sample_weights=test_sample_weights)}
        split_results.append({'scores': split_scores})

    return split_results


r_base = importr('base')
r_biobase = importr('Biobase')

metrics = ['roc_auc', 'average_precision', 'balanced_accuracy']
ordinal_encoder_categories = {
    'tumor_stage': ['NA', 'x', 'i', 'i or ii', 'ii', 'iii', 'iv']}

parser = ArgumentParser()
parser.add_argument('--data-dir', type=str, default='data', help='data dir')
parser.add_argument('--out-dir', type=str, default=os.getcwd(), help='out dir')
parser.add_argument('--n-jobs', type=int, default=-1, help='num parallel jobs')
parser.add_argument('--verbose', type=int, default=1, help='verbosity')
args = parser.parse_args()

all_X, all_y, all_groups, all_sample_weights, cv_split_flags = (
    [], [], [], [], [])
eset_files = sorted(glob('{}/tcga_*_resp_*_eset.rds'.format(args.data_dir))
                    + glob('{}/tcga_*_rest_*_eset.rds'.format(args.data_dir)))
num_esets = len(eset_files)
for eset_idx, eset_file in enumerate(eset_files):
    file_basename = os.path.splitext(os.path.split(eset_file)[1])[0]
    _, cancer, _, target, _, *rest = file_basename.split('_')
    cv_split_flag = cancer == 'stad' and target == 'oxaliplatin'

    if args.verbose < 2:
        print('Loading {:d}/{:d} esets'.format(eset_idx + 1, num_esets),
              end='\r', flush=True)
    else:
        print('Loading {}'.format(file_basename))

    eset = r_base.readRDS(eset_file)
    sample_meta = r_biobase.pData(eset)
    X = pd.DataFrame(index=sample_meta.index)
    y = np.array(sample_meta['Class'], dtype=int)

    if 'Group' in sample_meta.columns:
        groups = np.array(sample_meta['Group'], dtype=int)
        _, group_indices, group_counts = np.unique(
            groups, return_inverse=True, return_counts=True)
        sample_weights = (np.max(group_counts) / group_counts)[group_indices]
    else:
        groups = None
        sample_weights = None

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
    all_sample_weights.append(sample_weights)
    cv_split_flags.append(cv_split_flag)

if args.verbose < 2:
    print(flush=True)

print('Running SVM models')
all_results = Parallel(n_jobs=args.n_jobs, verbose=args.verbose)(
    delayed(fit_models)(X, y, groups, sample_weights, cv_split_flag)
    for X, y, groups, sample_weights, cv_split_flag in
    zip(all_X, all_y, all_groups, all_sample_weights, cv_split_flags))

if args.verbose < 1:
    print(flush=True)

mean_scores = []
all_roc_scores_df = None
all_pr_scores_df = None
for eset_file, split_results in zip(eset_files, all_results):
    file_basename = os.path.splitext(os.path.split(eset_file)[1])[0]
    _, cancer, analysis, target, data_type, *rest = file_basename.split('_')
    data_type = 'expr' if data_type == 'htseq' else data_type

    roc_scores, pr_scores = [], []
    for split_result in split_results:
        roc_scores.append(split_result['scores']['te']['roc_auc'])
        pr_scores.append(split_result['scores']['te']['pr_auc'])

    mean_score = np.mean(roc_scores)
    mean_scores.append([analysis, cancer, target, data_type, mean_score])

    dataset_name = '_'.join(file_basename.split('_')[:-1])
    model_name = '_'.join([dataset_name, 'svm'])
    dump(split_results, '{}/{}_split_results.pkl'.format(args.out_dir,
                                                         model_name))

    roc_scores_df = pd.DataFrame({dataset_name: roc_scores})
    pr_scores_df = pd.DataFrame({dataset_name: pr_scores})
    if all_roc_scores_df is None:
        all_roc_scores_df = roc_scores_df
        all_pr_scores_df = pr_scores_df
    else:
        all_roc_scores_df = pd.concat([all_roc_scores_df, roc_scores_df],
                                      axis=1)
        all_pr_scores_df = pd.concat([all_pr_scores_df, pr_scores_df], axis=1)

all_roc_scores_df.to_csv('{}/svm_covariate_model_roc_scores.tsv'
                         .format(args.out_dir), sep='\t')
all_pr_scores_df.to_csv('{}/svm_covariate_model_pr_scores.tsv'
                        .format(args.out_dir), sep='\t')

dump(all_roc_scores_df,
     '{}/svm_covariate_model_roc_scores.pkl'.format(args.out_dir))
dump(all_pr_scores_df,
     '{}/svm_covariate_model_pr_scores.pkl'.format(args.out_dir))

r_base.saveRDS(all_roc_scores_df,
               '{}/svm_covariate_model_roc_scores.rds'.format(args.out_dir))
r_base.saveRDS(all_pr_scores_df,
               '{}/svm_covariate_model_pr_scores.rds'.format(args.out_dir))

mean_scores_df = pd.DataFrame(mean_scores, columns=[
    'Analysis', 'Cancer', 'Target', 'Data Type', 'Mean Score'])
mean_scores_df.to_csv('{}/svm_covariate_model_mean_scores.tsv'
                      .format(args.out_dir), index=False, sep='\t')
if args.verbose > 0:
    print(tabulate(
        mean_scores_df.sort_values(
            ['Analysis', 'Cancer', 'Target', 'Data Type']),
        floatfmt='.4f', showindex=False, headers='keys'))