import os
import re
import warnings
from argparse import ArgumentParser

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
from joblib import dump, load
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr

numpy2ri.activate()
pandas2ri.activate()

parser = ArgumentParser()
parser.add_argument('--results-dir', type=str, default='results',
                    help='results dir')
args = parser.parse_args()

metric = {'surv': 'score', 'resp': 'roc_auc'}

r_base = importr('base')

all_scores_dfs = {}
split_results_regex = re.compile('^(.+?_(?:cnet|rfe))_split_results\\.pkl$')
for dirpath, dirnames, filenames in sorted(os.walk(args.results_dir)):
    for filename in filenames:
        if m := re.search(split_results_regex, filename):
            model_name = m.group(1)
            _, cancer, analysis, target, data_type, *rest = (
                model_name.split('_'))
            data_type = 'expr' if data_type == 'htseq' else data_type
            model_code = rest[-1]
            split_results_file = '{}/{}'.format(dirpath, filename)
            print('Loading', split_results_file)
            split_results = load(split_results_file)
            scores = []
            for split_result in split_results:
                if split_result is None:
                    scores.append(np.nan)
                else:
                    scores.append(split_result['scores']['te']
                                  [metric[analysis]])
            dataset_name = '_'.join(model_name.split('_')[:-1])
            scores_df = pd.DataFrame({dataset_name: scores})
            if model_code not in all_scores_dfs:
                all_scores_dfs[model_code] = scores_df
            else:
                all_scores_dfs[model_code] = pd.concat(
                    [all_scores_dfs[model_code], scores_df], axis=1)

for model_code, all_scores_df in all_scores_dfs.items():
    analysis = 'surv' if model_code == 'cnet' else 'resp'
    out_dir = '{}/{}'.format(args.results_dir, analysis)
    os.makedirs(out_dir, mode=0o755, exist_ok=True)
    all_scores_df.to_csv('{}/{}_model_scores.tsv'.format(out_dir, model_code),
                         sep='\t')
    dump(all_scores_df, '{}/{}_model_scores.pkl'.format(out_dir, model_code))
    r_base.saveRDS(all_scores_df, '{}/{}_model_scores.rds'
                   .format(out_dir, model_code))
