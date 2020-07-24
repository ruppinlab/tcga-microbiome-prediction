#!/bin/bash

export BLIS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

export MPLBACKEND=agg
if [ -d "/lscratch/$SLURM_JOB_ID" ]; then
    export TMPDIR=/lscratch/$SLURM_JOB_ID
else
    mkdir -p ~/tmp
    export TMPDIR=~/tmp
fi
export JOBLIB_TEMP_FOLDER=$TMPDIR
export PYTHONUNBUFFERED=1

python_warnings=(
    'ignore::FutureWarning:sklearn.utils.deprecation'
    'ignore::FutureWarning:rpy2.robjects.pandas2ri'
    'ignore:Optimization did not converge:UserWarning'
    'ignore:Optimization terminated early:UserWarning'
    'ignore:Persisting input arguments took:UserWarning'
    'ignore:Estimator fit failed:RuntimeWarning'
    'ignore:all coefficients are zero:UserWarning:sksurv_extensions.linear_model._coxnet'
    'ignore:divide by zero encountered in true_divide:RuntimeWarning:sksurv.linear_model.coxph'
    'ignore:invalid value encountered in true_divide:RuntimeWarning'
    'ignore:overflow encountered in exp:RuntimeWarning:sksurv.linear_model.coxph'
    'ignore:overflow encountered in power:RuntimeWarning:sksurv.linear_model.coxph'
)
OIFS="$IFS"
IFS=','
export PYTHONWARNINGS="${python_warnings[*]}"
IFS="$OIFS"

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate tcga-microbiome-prediction

if [[ -v SLURM_SUBMIT_DIR ]]; then
    SCRIPT_DIR="$SLURM_SUBMIT_DIR"
else
    SCRIPT_DIR=$(dirname $(realpath -s $0))
fi

RUN_MODEL_CMD="python $SCRIPT_DIR/run_model.py $@"
echo $RUN_MODEL_CMD
$RUN_MODEL_CMD
