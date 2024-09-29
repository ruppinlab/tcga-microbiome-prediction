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
export PYTHONUNBUFFERED=1

python_warnings=(
    'ignore:Optimization did not converge:UserWarning'
    'ignore:Optimization terminated early:UserWarning'
    'ignore:Persisting input arguments took:UserWarning'
    'ignore:Possible name collisions between functions:UserWarning'
    'ignore:A worker stopped while some jobs were given to the executor:UserWarning'
    'ignore:Estimator fit failed:RuntimeWarning'
    'ignore:Some fits failed:RuntimeWarning:sklearn_extensions.model_selection._validation'
    'ignore:Solver terminated early:UserWarning:sklearn.svm._base'
    'ignore:The max_iter was reached which means the coef_ did not converge:UserWarning:sklearn.linear_model._sag'
    'ignore:No features were selected:UserWarning:sklearn_extensions.feature_selection._base'
    'ignore:all coefficients are zero:UserWarning:sksurv_extensions.linear_model._coxnet'
    'ignore:divide by zero encountered in true_divide:RuntimeWarning:sksurv.linear_model.coxph'
    'ignore:divide by zero encountered in divide:RuntimeWarning:sksurv.linear_model.coxph'
    'ignore:invalid value encountered in true_divide:RuntimeWarning'
    'ignore:invalid value encountered in divide:RuntimeWarning'
    'ignore:overflow encountered in exp:RuntimeWarning:sksurv.linear_model.coxph'
    'ignore:overflow encountered in power:RuntimeWarning:sksurv.linear_model.coxph'
)
OIFS="$IFS"
IFS=','
export PYTHONWARNINGS="${python_warnings[*]}"
IFS="$OIFS"

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
while [[ -v CONDA_DEFAULT_ENV ]]; do
    conda deactivate
done
conda activate tcga-microbiome-prediction

if [[ -v SLURM_SUBMIT_DIR ]]; then
    SCRIPT_DIR="$SLURM_SUBMIT_DIR"
else
    SCRIPT_DIR=$(dirname $(realpath -s $0))
fi

RUN_MODEL_CMD="python $SCRIPT_DIR/run_model.py $@"
echo $RUN_MODEL_CMD
$RUN_MODEL_CMD
