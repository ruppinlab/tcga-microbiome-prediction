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

r_base = importr('base')

parser = ArgumentParser()
parser.add_argument('--model-code', type=str, choices=['cnet', 'rfe'],
                    help='model code')
parser.add_argument('--results-dir', type=str, default='results',
                    help='results dir')
parser.add_argument('--out-dir', type=str, default=os.getcwd(),
                    help='out dir')
args = parser.parse_args()

metric = {'surv': 'score', 'resp': 'roc_auc'}

all_scores_df = None
split_results_regex = re.compile('^(.+?)_split_results\\.pkl$')
for dirpath, dirnames, filenames in sorted(os.walk(args.results_dir)):
    for filename in filenames:
        if m := re.search(split_results_regex, filename):
            model_name = m.group(1)
            _, cancer, analysis, target, data_type, *rest = (
                model_name.split('_'))
            data_type = 'expr' if data_type == 'htseq' else data_type
            model_code = rest[-1]
            if model_code == args.model_code:
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
                if all_scores_df is None:
                    all_scores_df = scores_df
                else:
                    all_scores_df = pd.concat([all_scores_df, scores_df],
                                              axis=1)

all_scores_df.to_csv('{}/{}_model_scores.tsv'
                     .format(args.out_dir, args.model_code), sep='\t')
dump(all_scores_df, '{}/{}_model_scores.pkl'
     .format(args.out_dir, args.model_code))
r_base.saveRDS(all_scores_df, '{}/{}_model_scores.rds'
               .format(args.out_dir, args.model_code))
