import os
import sys
from argparse import ArgumentParser
from glob import glob
from shutil import rmtree

parser = ArgumentParser()
parser.add_argument('--data-dir', type=str, default='data')
args = parser.parse_args()

sbatch_opts = '--gres=lscratch:200 --mem-per-cpu=4096m --time=4-00:00:00'

project_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
eset_files = sorted(glob('{}/tcga_*_eset.rds'.format(args.data_dir)))
for eset_file in eset_files:
    dataset_name = os.path.splitext(os.path.split(eset_file)[1])[0]
    _, cancer, analysis, target, data_type, *rest = dataset_name.split('_')
    if analysis == 'surv':
        model_types = ['cnet']
    else:
        model_types = ['rfe', 'lgr']
        model_types.append('edger' if data_type == 'htseq' else 'limma')
    for model_type in model_types:
        model_name = '_'.join([dataset_name.rpartition('_')[0], model_type])
        results_dir = 'results/models/{}/{}'.format(analysis, model_name)
        if os.path.exists(results_dir):
            rmtree(results_dir)
        os.makedirs(results_dir, mode=0o755)
        num_jobs = 32 if data_type == 'kraken' else 56
        cmd = ('{}/submit_slurm_model.sh --model-type {} --dataset {} '
               '--n-jobs {:d} --sbatch-opts "{} --output={}/slurm-\\%j.out"'
               .format(project_dir, model_type, eset_file, num_jobs,
                       sbatch_opts, results_dir))
        print(cmd)
        os.system(cmd)
