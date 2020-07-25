import os
import sys
from argparse import ArgumentParser
from glob import glob
from shutil import rmtree

parser = ArgumentParser()
parser.add_argument('--data-dir', type=str, default='data')
args = parser.parse_args()

sbatch_opts = '--gres=lscratch:400 --mem-per-cpu=3584m --time=4-00:00:00'

project_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
eset_files = sorted(glob('{}/tcga_*_eset.rds'.format(args.data_dir)))
for eset_file in eset_files:
    file_basename = os.path.splitext(os.path.split(eset_file)[1])[0]
    file_basename_parts = file_basename.split('_')
    data_type = file_basename_parts[4]
    num_jobs = 32 if data_type == 'kraken' else 56
    results_dir = 'results/{}/{}'.format(
        file_basename_parts[2],
        '_'.join(file_basename_parts[:-1] + ['cnet']))
    if os.path.exists(results_dir):
        rmtree(results_dir)
    os.makedirs(results_dir, mode=0o755)
    cmd = ('{}/submit_slurm_model.sh --dataset {} --n-jobs {:d} '
           '--sbatch-opts "{} --output={}/slurm-\\%j.out"'
           .format(project_dir, eset_file, num_jobs, sbatch_opts, results_dir))
    print(cmd)
    os.system(cmd)
