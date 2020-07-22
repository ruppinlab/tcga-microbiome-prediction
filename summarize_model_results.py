#!/usr/bin/env python

import os
import re
from argparse import ArgumentParser
from glob import glob

import numpy as np
import pandas as pd
from joblib import load

parser = ArgumentParser()
parser.add_argument('--cox-covar-file', type=str,
                    default='cox/cox_covariate_model_mean_scores.tsv',
                    help='cox covariate mean scores file')
parser.add_argument('--svm-covar-file', type=str,
                    default='svm/svm_covariate_model_mean_scores.tsv',
                    help='svm covariate mean scores file')
parser.add_argument('--results-dir', type=str, default='results',
                    help='results dir')
parser.add_argument('--out-dir', type=str, default=os.getcwd(), help='out dir')
args = parser.parse_args()

metric = {'surv': 'score', 'resp': 'roc_auc', 'rest': 'roc_auc'}
penalty_factor_meta_col = 'Penalty Factor'

results = []
split_results_regex = re.compile('^(.+?)_split_results\\.pkl$')
for dirpath, dirnames, filenames in sorted(os.walk(args.results_dir)):
    for filename in filenames:
        if m := re.search(split_results_regex, filename):
            model_name = m.group(1)
            _, cancer, analysis, target, data_type, *rest = (
                model_name.split('_'))
            data_type = 'expr' if data_type == 'htseq' else data_type
            model_code = rest[-1]

            slurm_file = glob('{}/slurm-*.out'.format(dirpath))[0]
            slurm_basename = os.path.splitext(os.path.split(slurm_file)[1])[0]
            job_id = re.findall('\\d+', slurm_basename)[0]

            split_results_file = '{}/{}'.format(dirpath, filename)
            print('Loading', split_results_file)
            split_results = load(split_results_file)

            scores = []
            num_features = []
            for split_result in split_results:
                if split_result is None:
                    continue
                scores.append(split_result['scores']['te'][metric[analysis]])
                split_feature_meta = split_result['feature_meta']
                if penalty_factor_meta_col in split_feature_meta.columns:
                    num_features.append(split_feature_meta.loc[
                        split_feature_meta[penalty_factor_meta_col]
                        != 0].shape[0])
                else:
                    num_features.append(split_feature_meta.shape[0])

            results.append([analysis, cancer, target, data_type, model_code,
                            np.mean(scores), np.mean(num_features), job_id])

results_summary = pd.DataFrame(results, columns=[
    'Analysis', 'Cancer', 'Target', 'Data Type', 'Model Code',
    'Mean Score', 'Mean Num Features', 'Job ID'])
if os.path.isfile(args.cox_covar_file):
    cox_mean_scores = pd.read_csv(args.cox_covar_file, sep='\t')
    results_summary = pd.merge(
        results_summary, cox_mean_scores, how='left',
        on=['Analysis', 'Cancer', 'Target', 'Data Type'])
if os.path.isfile(args.svm_covar_file):
    svm_mean_scores = pd.read_csv(args.svm_covar_file, sep='\t')
    results_summary = pd.merge(
        results_summary, svm_mean_scores, how='left',
        on=['Analysis', 'Cancer', 'Target', 'Data Type'])
results_summary['Covariate Mean Score'] = (
    results_summary['Mean Score_y'].combine_first(
        results_summary['Mean Score']))
results_summary.drop(columns=['Mean Score', 'Mean Score_y'], inplace=True)
results_summary.rename(columns={'Mean Score_x': 'Mean Score'}, inplace=True)
results_summary.sort_values(by=['Analysis', 'Cancer', 'Target', 'Data Type',
                                'Model Code', 'Mean Score'], inplace=True)
results_summary.to_csv('{}/model_results_summary.tsv'.format(args.out_dir),
                       sep='\t', index=False)
