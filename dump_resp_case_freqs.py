import os
import warnings
from argparse import ArgumentParser
from glob import glob

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
from tabulate import tabulate

numpy2ri.activate()
pandas2ri.activate()

parser = ArgumentParser()
parser.add_argument('--data-dir', type=str, default='data', help='data dir')
parser.add_argument('--results-dir', type=str, default='results/models',
                    help='results dir')
parser.add_argument('--sort-by', type=str, nargs='+',
                    choices=['Dataset Name', 'Analysis', 'Num Cases',
                             'NR (-) Cases', 'R (+) Cases'],
                    default=['Analysis', 'Dataset Name'],
                    help='columns to sort table by')
args = parser.parse_args()

analysis = 'resp'

r_base = importr('base')
r_biobase = importr('Biobase')

results = []
eset_files = glob('{}/tcga_*_resp_*_eset.rds'.format(args.data_dir))
for eset_file in sorted(eset_files):
    file_basename = os.path.splitext(os.path.split(eset_file)[1])[0]
    _, cancer, analysis, target, data_type, *rest = file_basename.split('_')
    dataset_name = '_'.join(file_basename.split('_')[:-1])
    eset = r_base.readRDS(eset_file)
    sample_meta = r_biobase.pData(eset)
    num_cases = sample_meta['case_submitter_id'].nunique()
    neg_cases, pos_cases = (sample_meta.groupby('Class')['case_submitter_id']
                            .nunique())
    results.append([dataset_name, analysis, num_cases, neg_cases, pos_cases])

results_df = pd.DataFrame(results, columns=[
    'Dataset Name', 'Analysis', 'Num Cases', 'NR (-) Cases', 'R (+) Cases'])
out_dir = '{}/{}'.format(args.results_dir, analysis)
os.makedirs(out_dir, mode=0o755, exist_ok=True)
results_df.to_csv('{}/resp_case_freqs.tsv'.format(out_dir), sep='\t',
                  index=False)
print(tabulate(results_df.sort_values(by=args.sort_by), headers='keys',
               showindex=False))
