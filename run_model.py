import atexit
import os
import re
import sys
import warnings
from argparse import ArgumentParser
from decimal import Decimal
from glob import glob
from pprint import pprint
from shutil import rmtree
from tempfile import mkdtemp, gettempdir
from traceback import format_exception_only

warnings.filterwarnings('ignore', category=FutureWarning,
                        module='sklearn.utils.deprecation')
warnings.filterwarnings('ignore', category=FutureWarning,
                        module='rpy2.robjects.pandas2ri')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rpy2.rinterface_lib.embedded as r_embedded

r_embedded.set_initoptions(
    ('rpy2', '--quiet', '--no-save', '--max-ppsize=500000'))

import rpy2.robjects as robjects
import seaborn as sns
from eli5 import explain_weights_df
from joblib import Memory, Parallel, delayed, dump, parallel_backend
from matplotlib.offsetbox import AnchoredText
from pandas.api.types import (
    is_bool_dtype, is_categorical_dtype, is_integer_dtype, is_float_dtype,
    is_object_dtype, is_string_dtype)
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr
from scipy.stats import iqr
from sklearn.base import BaseEstimator, clone
from sklearn.compose import ColumnTransformer
from sklearn.exceptions import FitFailedWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (auc, average_precision_score,
                             balanced_accuracy_score, precision_recall_curve,
                             roc_auc_score, roc_curve)
from sklearn.model_selection import ParameterGrid, RepeatedStratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import (OneHotEncoder, OrdinalEncoder,
                                   StandardScaler)
from sklearn.svm import SVC
from sklearn.utils import check_random_state, _determine_key_type
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv
from tabulate import tabulate

numpy2ri.activate()
pandas2ri.activate()

from sklearn_extensions.compose import ExtendedColumnTransformer
from sklearn_extensions.feature_selection import (
    EdgeR, EdgeRFilterByExpr, ExtendedRFE, Limma, SelectFromModel)
from sklearn_extensions.linear_model import CachedLogisticRegression
from sklearn_extensions.model_selection import (
    ExtendedGridSearchCV, RepeatedStratifiedGroupKFold, shuffle_y)
from sklearn_extensions.pipeline import (ExtendedPipeline,
                                         transform_feature_meta)
from sklearn_extensions.preprocessing import EdgeRTMMLogCPM
from sksurv_extensions.model_selection import (
    RepeatedSurvivalStratifiedKFold, SurvivalStratifiedShuffleSplit,
    SurvivalStratifiedSampleFromGroupShuffleSplit)
from sksurv_extensions.linear_model import (
    CachedExtendedCoxnetSurvivalAnalysis, MetaCoxnetSurvivalAnalysis)


def warning_format(message, category, filename, lineno, file=None, line=None):
    return ' {}: {}'.format(category.__name__, message)


def load_dataset(dataset_file):
    dataset_name, file_extension = os.path.splitext(
        os.path.split(dataset_file)[1])
    if not os.path.isfile(dataset_file) or file_extension.lower() != '.rds':
        raise IOError('File does not exist/invalid: {}'.format(dataset_file))
    eset = r_base.readRDS(dataset_file)
    X = pd.DataFrame(r_base.t(r_biobase.exprs(eset)),
                     columns=r_biobase.featureNames(eset),
                     index=r_biobase.sampleNames(eset))
    sample_meta = r_biobase.pData(eset)
    if analysis == 'surv':
        y = Surv.from_dataframe('Status', 'Survival_in_days', sample_meta)
    else:
        y = np.array(sample_meta['Class'], dtype=int)
    if 'Group' in sample_meta.columns:
        groups = np.array(sample_meta['Group'], dtype=int)
        _, group_indices, group_counts = np.unique(
            groups, return_inverse=True, return_counts=True)
        if ('GroupWeight' in sample_meta.columns
                and sample_meta['GroupWeight'].unique().size > 1):
            group_weights = np.array(sample_meta['GroupWeight'], dtype=float)
        else:
            group_weights = None
        sample_weights = (np.max(group_counts) / group_counts)[group_indices]
    else:
        groups = None
        group_weights = None
        sample_weights = None
    try:
        feature_meta = r_biobase.fData(eset)
        feature_meta_category_cols = (
            feature_meta.select_dtypes(include='category').columns)
        feature_meta[feature_meta_category_cols] = (
            feature_meta[feature_meta_category_cols].astype(str))
    except ValueError:
        feature_meta = pd.DataFrame(index=r_biobase.featureNames(eset))
    new_feature_names = []
    if penalty_factor_meta_col in feature_meta.columns:
        raise RuntimeError('{} column already exists in feature_meta'
                           .format(penalty_factor_meta_col))
    feature_meta[penalty_factor_meta_col] = 1
    for sample_meta_col in sample_meta_cols:
        if sample_meta_col not in sample_meta.columns:
            raise RuntimeError('{} column does not exist in sample_meta'
                               .format(sample_meta_col))
        if sample_meta_col in X.columns:
            raise RuntimeError('{} column already exists in X'
                               .format(sample_meta_col))
        is_category = (is_categorical_dtype(sample_meta[sample_meta_col])
                       or is_object_dtype(sample_meta[sample_meta_col])
                       or is_string_dtype(sample_meta[sample_meta_col]))
        if not is_category:
            X[sample_meta_col] = sample_meta[sample_meta_col]
            new_feature_names.append(sample_meta_col)
        elif sample_meta_col in ordinal_encoder_categories:
            if sample_meta_col not in ordinal_encoder_categories:
                raise RuntimeError('No ordinal encoder categories config '
                                   'exists for {}'.format(sample_meta_col))
            if sample_meta[sample_meta_col].unique().size > 1:
                ode = OrdinalEncoder(categories=[
                    ordinal_encoder_categories[sample_meta_col]])
                ode.fit(sample_meta[[sample_meta_col]])
                X[sample_meta_col] = ode.transform(
                    sample_meta[[sample_meta_col]])
                new_feature_names.append(sample_meta_col)
        else:
            num_categories = sample_meta[sample_meta_col][
                sample_meta[sample_meta_col] != 'NA'].unique().size
            if num_categories > 2:
                ohe_drop = (['NA'] if 'NA' in
                            sample_meta[sample_meta_col].values else None)
                ohe = OneHotEncoder(drop=ohe_drop, sparse=False)
                ohe.fit(sample_meta[[sample_meta_col]])
                new_sample_meta_cols = []
                for category in ohe.categories_[0]:
                    if category == 'NA':
                        continue
                    new_sample_meta_col = '{}_{}'.format(
                        sample_meta_col, category).replace(' ', '_')
                    new_sample_meta_cols.append(new_sample_meta_col)
                X = X.join(pd.DataFrame(
                    ohe.transform(sample_meta[[sample_meta_col]]),
                    index=sample_meta[[sample_meta_col]].index,
                    columns=new_sample_meta_cols), sort=False)
                new_feature_names.extend(new_sample_meta_cols)
            elif num_categories == 2:
                ohe = OneHotEncoder(drop='first', sparse=False)
                ohe.fit(sample_meta[[sample_meta_col]])
                category = ohe.categories_[0][1]
                new_sample_meta_col = '{}_{}'.format(
                    sample_meta_col, category).replace(' ', '_')
                X[new_sample_meta_col] = ohe.transform(
                    sample_meta[[sample_meta_col]])
                new_feature_names.append(new_sample_meta_col)
    new_feature_meta = pd.DataFrame(index=new_feature_names)
    for feature_meta_col in feature_meta.columns:
        if (is_categorical_dtype(feature_meta[feature_meta_col])
                or is_object_dtype(feature_meta[feature_meta_col])
                or is_string_dtype(feature_meta[feature_meta_col])):
            new_feature_meta[feature_meta_col] = ''
        elif (is_integer_dtype(feature_meta[feature_meta_col])
              or is_float_dtype(feature_meta[feature_meta_col])):
            new_feature_meta[feature_meta_col] = 0
        elif is_bool_dtype(feature_meta[feature_meta_col]):
            new_feature_meta[feature_meta_col] = False
    new_feature_meta[penalty_factor_meta_col] = 0
    feature_meta = feature_meta.append(new_feature_meta, verify_integrity=True)
    return (dataset_name, X, y, groups, group_weights, sample_weights,
            sample_meta, feature_meta)


def get_col_trf_col_grps(X, col_trf_pat_grps):
    X_ct = X.copy()
    col_trf_col_grps = []
    for col_trf_pats in col_trf_pat_grps:
        col_trf_cols = []
        for pattern in col_trf_pats:
            col_trf_cols.append(X_ct.columns.str.contains(pattern, regex=True))
        X_ct = X_ct.loc[:, col_trf_cols[0]]
        col_trf_col_grps.append(col_trf_cols)
    return col_trf_col_grps


def setup_pipe_and_param_grid(X):
    clf_c = np.logspace(-5, 3, 9)
    l1_ratio = np.array([0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 1.])
    skb_k = np.insert(np.linspace(2, 400, num=200, dtype=int), 0, 1)
    sfm_c = (np.logspace(-2, 3, 6) if data_type == 'kraken' else
             np.logspace(-2, 1, 4))
    # prognosis: coxnet
    if analysis == 'surv':
        if data_type == 'kraken':
            pipe = ExtendedPipeline(
                memory=memory,
                param_routing={'srv1': ['feature_meta']},
                steps=[
                    ('trf0', StandardScaler()),
                    ('srv1', MetaCoxnetSurvivalAnalysis(
                        estimator=CachedExtendedCoxnetSurvivalAnalysis(
                            alpha_min_ratio=0.01, fit_baseline_model=True,
                            max_iter=1000000, memory=memory, n_alphas=100,
                            penalty_factor_meta_col=penalty_factor_meta_col,
                            normalize=False, penalty_factor=None)))])
            param_grid_dict = {'srv1__estimator__l1_ratio': l1_ratio}
        else:
            col_trf_col_grps = get_col_trf_col_grps(X, [['^ENSG.+$']])
            pipe = ExtendedPipeline(
                memory=memory,
                param_routing={'srv2': ['feature_meta'],
                               'trf0': ['sample_meta']},
                steps=[
                    ('trf0', ExtendedColumnTransformer(
                        n_jobs=1,
                        param_routing={'trf0': ['sample_meta']},
                        remainder='passthrough',
                        transformers=[
                            ('trf0', ExtendedPipeline(
                                memory=memory,
                                param_routing={'slr0': ['sample_meta'],
                                               'trf1': ['sample_meta']},
                                steps=[
                                    ('slr0', EdgeRFilterByExpr(
                                        is_classif=False)),
                                    ('trf1', EdgeRTMMLogCPM(
                                        prior_count=1))]),
                             col_trf_col_grps[0][0])])),
                    ('trf1', StandardScaler()),
                    ('srv2', MetaCoxnetSurvivalAnalysis(
                        estimator=CachedExtendedCoxnetSurvivalAnalysis(
                            alpha_min_ratio=0.01, fit_baseline_model=True,
                            max_iter=1000000, memory=memory, n_alphas=100,
                            penalty_factor_meta_col=penalty_factor_meta_col,
                            normalize=False, penalty_factor=None)))])
            param_grid_dict = {'srv2__estimator__l1_ratio': l1_ratio}
    # drug response: svm-rfe
    elif args.model_type == 'rfe':
        if data_type == 'kraken':
            pipe = ExtendedPipeline(
                memory=memory,
                param_routing={'clf1': ['feature_meta', 'sample_weight']},
                steps=[
                    ('trf0', StandardScaler()),
                    ('clf1', ExtendedRFE(
                        estimator=SVC(
                            cache_size=2000, class_weight='balanced',
                            kernel='linear', max_iter=int(1e8),
                            random_state=random_seed),
                        memory=memory, n_features_to_select=None,
                        penalty_factor_meta_col=penalty_factor_meta_col,
                        reducing_step=False, step=1, tune_step_at=None,
                        tuning_step=1))])
            param_grid_dict = {'clf1__estimator__C': clf_c,
                               'clf1__n_features_to_select': skb_k}
        else:
            col_trf_col_grps = get_col_trf_col_grps(X, [['^ENSG.+$']])
            pipe = ExtendedPipeline(
                memory=memory,
                param_routing={'clf2': ['feature_meta', 'sample_weight'],
                               'trf0': ['sample_meta']},
                steps=[
                    ('trf0', ExtendedColumnTransformer(
                        n_jobs=1,
                        param_routing={'trf0': ['sample_meta']},
                        remainder='passthrough',
                        transformers=[
                            ('trf0', ExtendedPipeline(
                                memory=memory,
                                param_routing={'slr0': ['sample_meta'],
                                               'trf1': ['sample_meta']},
                                steps=[
                                    ('slr0', EdgeRFilterByExpr(
                                        is_classif=True)),
                                    ('trf1', EdgeRTMMLogCPM(
                                        prior_count=1))]),
                             col_trf_col_grps[0][0])])),
                    ('trf1', StandardScaler()),
                    ('clf2', ExtendedRFE(
                        estimator=SVC(
                            cache_size=2000, class_weight='balanced',
                            kernel='linear', max_iter=int(1e8),
                            random_state=random_seed),
                        memory=memory, n_features_to_select=None,
                        penalty_factor_meta_col=penalty_factor_meta_col,
                        reducing_step=True, step=0.05, tune_step_at=1300,
                        tuning_step=1))])
            param_grid_dict = {'clf2__estimator__C': clf_c,
                               'clf2__n_features_to_select': skb_k}
    # drug response: elasticnet logistic regression
    elif args.model_type == 'lgr':
        if data_type == 'kraken':
            col_trf_col_grps = get_col_trf_col_grps(
                X, [['^(?!gender_male|age_at_diagnosis|tumor_stage).*$',
                     '^(?:gender_male|age_at_diagnosis|tumor_stage)$']])
            pipe = ExtendedPipeline(
                memory=memory,
                param_routing={'clf1': ['sample_weight'],
                               'trf0': ['sample_weight']},
                steps=[
                    ('trf0', ExtendedColumnTransformer(
                        n_jobs=1,
                        param_routing={'trf0': ['sample_weight']},
                        remainder='passthrough',
                        transformers=[
                            ('trf0', ExtendedPipeline(
                                memory=memory,
                                param_routing={'slr1': ['sample_weight']},
                                steps=[
                                    ('trf0', StandardScaler()),
                                    ('slr1', SelectFromModel(
                                        estimator=CachedLogisticRegression(
                                            class_weight='balanced',
                                            max_iter=5000,
                                            memory=memory,
                                            penalty='elasticnet',
                                            random_state=random_seed,
                                            solver='saga'),
                                        max_features=400,
                                        threshold=1e-10))]),
                             col_trf_col_grps[0][0]),
                            ('trf1', ExtendedPipeline(
                                memory=memory,
                                param_routing=None,
                                steps=[('trf0', StandardScaler())]),
                             col_trf_col_grps[0][1])])),
                    ('clf1', LogisticRegression(
                        class_weight='balanced',
                        max_iter=5000,
                        penalty='l2',
                        random_state=random_seed,
                        solver='saga'))])
            param_grid_dict = {
                'clf1__C': clf_c,
                'trf0__trf0__slr1__estimator__C': sfm_c,
                'trf0__trf0__slr1__estimator__l1_ratio': l1_ratio}
        elif data_type == 'htseq':
            col_trf_col_grps = get_col_trf_col_grps(
                X, [['^ENSG.+$', '^(?!ENSG).*$']])
            pipe = ExtendedPipeline(
                memory=memory,
                param_routing={'clf1': ['sample_weight'],
                               'trf0': ['sample_meta', 'sample_weight']},
                steps=[
                    ('trf0', ExtendedColumnTransformer(
                        n_jobs=1,
                        param_routing={'trf0': ['sample_meta',
                                                'sample_weight']},
                        remainder='passthrough',
                        transformers=[
                            ('trf0', ExtendedPipeline(
                                memory=memory,
                                param_routing={'slr0': ['sample_meta'],
                                               'slr3': ['sample_weight'],
                                               'trf1': ['sample_meta']},
                                steps=[
                                    ('slr0', EdgeRFilterByExpr(
                                        is_classif=True)),
                                    ('trf1', EdgeRTMMLogCPM(
                                        prior_count=1)),
                                    ('trf2', StandardScaler()),
                                    ('slr3', SelectFromModel(
                                        estimator=CachedLogisticRegression(
                                            class_weight='balanced',
                                            max_iter=5000,
                                            memory=memory,
                                            penalty='elasticnet',
                                            random_state=random_seed,
                                            solver='saga'),
                                        max_features=400,
                                        threshold=1e-10))]),
                             col_trf_col_grps[0][0]),
                            ('trf1', ExtendedPipeline(
                                memory=memory,
                                param_routing=None,
                                steps=[('trf0', StandardScaler())]),
                             col_trf_col_grps[0][1])])),
                    ('clf1', LogisticRegression(
                        class_weight='balanced',
                        max_iter=5000,
                        penalty='l2',
                        random_state=random_seed,
                        solver='saga'))])
            param_grid_dict = {
                'clf1__C': clf_c,
                'trf0__trf0__slr3__estimator__C': sfm_c,
                'trf0__trf0__slr3__estimator__l1_ratio': l1_ratio}
        else:
            col_trf_col_grps = get_col_trf_col_grps(
                X, [['^(?!gender_male|age_at_diagnosis|tumor_stage).*$',
                     '^(?:gender_male|age_at_diagnosis|tumor_stage)$'],
                    ['^ENSG.+$']])
            pipe = ExtendedPipeline(
                memory=memory,
                param_routing={'clf1': ['sample_weight'],
                               'trf0': ['sample_meta', 'sample_weight']},
                steps=[
                    ('trf0', ExtendedColumnTransformer(
                        n_jobs=1,
                        param_routing={'trf0': ['sample_meta',
                                                'sample_weight']},
                        remainder='passthrough',
                        transformers=[
                            ('trf0', ExtendedPipeline(
                                memory=memory,
                                param_routing={'slr2': ['sample_weight'],
                                               'trf0': ['sample_meta']},
                                steps=[
                                    ('trf0', ExtendedColumnTransformer(
                                        n_jobs=1,
                                        param_routing={
                                            'trf0': ['sample_meta']},
                                        remainder='passthrough',
                                        transformers=[
                                            ('trf0', ExtendedPipeline(
                                                memory=memory,
                                                param_routing={
                                                    'slr0': ['sample_meta'],
                                                    'trf1': ['sample_meta']},
                                                steps=[
                                                    ('slr0', EdgeRFilterByExpr(
                                                        is_classif=True)),
                                                    ('trf1', EdgeRTMMLogCPM(
                                                        prior_count=1))]),
                                             col_trf_col_grps[1][0])])),
                                    ('trf1', StandardScaler()),
                                    ('slr2', SelectFromModel(
                                        estimator=CachedLogisticRegression(
                                            class_weight='balanced',
                                            max_iter=5000,
                                            memory=memory,
                                            penalty='elasticnet',
                                            random_state=random_seed,
                                            solver='saga'),
                                        max_features=400,
                                        threshold=1e-10))]),
                             col_trf_col_grps[0][0]),
                            ('trf1', ExtendedPipeline(
                                memory=memory,
                                param_routing=None,
                                steps=[('trf0', StandardScaler())]),
                             col_trf_col_grps[0][1])])),
                    ('clf1', LogisticRegression(
                        class_weight='balanced',
                        max_iter=5000,
                        penalty='l2',
                        random_state=random_seed,
                        solver='saga'))])
            param_grid_dict = {
                'clf1__C': clf_c,
                'trf0__trf0__slr2__estimator__C': sfm_c,
                'trf0__trf0__slr2__estimator__l1_ratio': l1_ratio}
    # drug response: limma/edgeR L2 logistic regression
    elif data_type == 'kraken':
        col_trf_col_grps = get_col_trf_col_grps(
            X, [['^(?!gender_male|age_at_diagnosis|tumor_stage).*$']])
        pipe = ExtendedPipeline(
            memory=memory,
            param_routing={'clf2': ['sample_weight'],
                           'trf0': ['sample_meta']},
            steps=[
                ('trf0', ExtendedColumnTransformer(
                    n_jobs=1,
                    param_routing={'trf0': ['sample_meta']},
                    remainder='passthrough',
                    transformers=[
                        ('trf0', ExtendedPipeline(
                            memory=memory,
                            param_routing={'slr0': ['sample_meta']},
                            steps=[
                                ('slr0', Limma(
                                    memory=memory,
                                    robust=True,
                                    trend=True))]),
                         col_trf_col_grps[0][0])])),
                ('trf1', StandardScaler()),
                ('clf2', LogisticRegression(
                    class_weight='balanced',
                    max_iter=5000,
                    penalty='l2',
                    random_state=random_seed,
                    solver='saga'))])
        param_grid_dict = {'clf2__C': clf_c,
                           'trf0__trf0__slr0__k': skb_k}
    elif data_type == 'htseq':
        col_trf_col_grps = get_col_trf_col_grps(X, [['^ENSG.+$']])
        pipe = ExtendedPipeline(
            memory=memory,
            param_routing={'clf2': ['sample_weight'],
                           'trf0': ['sample_meta']},
            steps=[
                ('trf0', ExtendedColumnTransformer(
                    n_jobs=1,
                    param_routing={'trf0': ['sample_meta']},
                    remainder='passthrough',
                    transformers=[
                        ('trf0', ExtendedPipeline(
                            memory=memory,
                            param_routing={'slr0': ['sample_meta'],
                                           'slr1': ['sample_meta']},
                            steps=[
                                ('slr0', EdgeRFilterByExpr(
                                    is_classif=True)),
                                ('slr1', EdgeR(
                                    memory=memory,
                                    prior_count=1,
                                    robust=True))]),
                         col_trf_col_grps[0][0])])),
                ('trf1', StandardScaler()),
                ('clf2', LogisticRegression(
                    class_weight='balanced',
                    max_iter=5000,
                    penalty='l2',
                    random_state=random_seed,
                    solver='saga'))])
        param_grid_dict = {'clf2__C': clf_c,
                           'trf0__trf0__slr1__k': skb_k}
    else:
        col_trf_col_grps = get_col_trf_col_grps(
            X, [['^(?!gender_male|age_at_diagnosis|tumor_stage).*$'],
                ['^ENSG.+$']])
        pipe = ExtendedPipeline(
            memory=memory,
            param_routing={'clf2': ['sample_weight'],
                           'trf0': ['sample_meta']},
            steps=[('trf0', ExtendedColumnTransformer(
                n_jobs=1,
                param_routing={'trf0': ['sample_meta']},
                remainder='passthrough',
                transformers=[
                    ('trf0', ExtendedPipeline(
                        memory=memory,
                        param_routing={'slr1': ['sample_meta'],
                                       'trf0': ['sample_meta']},
                        steps=[
                            ('trf0', ExtendedColumnTransformer(
                                n_jobs=1,
                                param_routing={'trf0': ['sample_meta']},
                                remainder='passthrough',
                                transformers=[
                                    ('trf0', ExtendedPipeline(
                                        memory=memory,
                                        param_routing={
                                            'slr0': ['sample_meta'],
                                            'trf1': ['sample_meta']},
                                        steps=[
                                            ('slr0', EdgeRFilterByExpr(
                                                is_classif=True)),
                                            ('trf1', EdgeRTMMLogCPM(
                                                prior_count=1))]),
                                     col_trf_col_grps[1][0])])),
                            ('slr1', Limma(
                                memory=memory,
                                robust=True,
                                trend=True))]),
                     col_trf_col_grps[0][0])])),
                   ('trf1', StandardScaler()),
                   ('clf2', LogisticRegression(
                       class_weight='balanced',
                       max_iter=5000,
                       penalty='l2',
                       random_state=random_seed,
                       solver='saga'))])
        param_grid_dict = {'clf2__C': clf_c,
                           'trf0__trf0__slr1__k': skb_k}
    param_grid = [param_grid_dict.copy()]
    return pipe, param_grid, param_grid_dict


def col_trf_info(col_trf):
    col_trf_col_strs = []
    for trf_name, trf_transformer, trf_cols in col_trf.transformers:
        col_trf_col_strs.append('{}: {:d}'.format(
            trf_name, (np.count_nonzero(trf_cols)
                       if _determine_key_type(trf_cols) == 'bool'
                       else trf_cols.shape[0])))
        if (isinstance(trf_transformer, Pipeline)
                and isinstance(trf_transformer[0], ColumnTransformer)):
            col_trf_col_strs.append(col_trf_info(trf_transformer[0]))
    return '({})'.format(' '.join(col_trf_col_strs))


def get_param_type(param):
    pipe_step_type_regex = re.compile(
        r'^({})\d+$'.format('|'.join(pipe_step_types)))
    param_parts = param.split('__')
    param_parts_start_idx = [i for i, p in enumerate(param_parts)
                             if pipe_step_type_regex.match(p)][-1]
    param_parts[param_parts_start_idx] = pipe_step_type_regex.sub(
        r'\1', param_parts[param_parts_start_idx])
    param_type = '__'.join(param_parts[param_parts_start_idx:])
    return param_type


def fit_pipeline(X, y, steps, params=None, param_routing=None,
                 fit_params=None):
    pipe = ExtendedPipeline(steps, memory=memory, param_routing=param_routing)
    if params is None:
        params = {}
    pipe.set_params(**params)
    if fit_params is None:
        fit_params = {}
    try:
        pipe.fit(X, y, **fit_params)
    except ArithmeticError as e:
        warnings.formatwarning = warning_format
        warnings.warn('Estimator fit failed. Details: {}'
                      .format(format_exception_only(type(e), e)[0]),
                      category=FitFailedWarning)
        pipe = None
    if args.scv_verbose == 0:
        print('.' if pipe is not None else 'x', end='', flush=True)
    return pipe


def calculate_test_scores(estimator, X_test, y_test, metrics,
                          predict_params=None, score_params=None):
    scores = {}
    if predict_params is None:
        predict_params = {}
    if hasattr(estimator, 'decision_function'):
        y_score = estimator.decision_function(X_test, **predict_params)
        scores['y_score'] = y_score
    y_pred = estimator.predict(X_test, **predict_params)
    scores['y_pred'] = y_pred
    if score_params is None:
        score_params = {}
    if isinstance(metrics, str):
        metrics = [metrics]
    for metric in metrics:
        if metric in ('concordance_index_censored', 'score'):
            scores[metric] = concordance_index_censored(
                y_test[y_test.dtype.names[0]], y_test[y_test.dtype.names[1]],
                y_pred)[0]
        elif metric == 'roc_auc':
            scores[metric] = roc_auc_score(
                y_test, y_score, **score_params)
            scores['fpr'], scores['tpr'], _ = roc_curve(
                y_test, y_score, pos_label=1, **score_params)
        elif metric == 'balanced_accuracy':
            scores[metric] = balanced_accuracy_score(
                y_test, y_pred, **score_params)
        elif metric == 'average_precision':
            scores[metric] = average_precision_score(
                y_test, y_score, **score_params)
            scores['pre'], scores['rec'], _ = precision_recall_curve(
                y_test, y_score, pos_label=1, **score_params)
            scores['pr_auc'] = auc(scores['rec'], scores['pre'])
    return scores


def get_perm_test_split_data(X, perm_y, cv, cv_params=None):
    if cv_params is None:
        cv_params = {}
    perm_split_idxs = list(cv.split(X, perm_y, **cv_params))
    return perm_y, perm_split_idxs


def fit_and_score(estimator, X_train, y_train, X_test, y_test, scoring,
                  fit_params=None, predict_params=None, score_params=None):
    if fit_params is None:
        fit_params = {}
    estimator.fit(X_train, y_train, **fit_params)
    scores = calculate_test_scores(estimator, X_test, y_test, scoring,
                                   predict_params=predict_params,
                                   score_params=score_params)
    return scores[scoring]


def get_final_feature_meta(pipe, feature_meta):
    for estimator in pipe:
        feature_meta = transform_feature_meta(estimator, feature_meta)
    final_estimator = pipe[-1]
    if isinstance(final_estimator, MetaCoxnetSurvivalAnalysis):
        feature_weights = final_estimator.coef_
        feature_weights = np.ravel(feature_weights)
        feature_mask = feature_weights != 0
        if penalty_factor_meta_col in feature_meta.columns:
            feature_mask[feature_meta[penalty_factor_meta_col] == 0] = True
        feature_meta = feature_meta.copy()
        feature_meta = feature_meta.loc[feature_mask]
        feature_meta['Weight'] = feature_weights[feature_mask]
    else:
        feature_weights = explain_weights_df(
            final_estimator, feature_names=feature_meta.index.values)
        if feature_weights is None and hasattr(final_estimator, 'estimator_'):
            feature_weights = explain_weights_df(
                final_estimator.estimator_,
                feature_names=feature_meta.index.values)
        if feature_weights is not None:
            feature_weights.set_index('feature', inplace=True,
                                      verify_integrity=True)
            feature_weights.columns = map(str.title, feature_weights.columns)
            feature_meta = feature_meta.join(feature_weights, how='inner')
            if (feature_meta['Weight'] == 0).any():
                if penalty_factor_meta_col in feature_meta.columns:
                    feature_meta = feature_meta.loc[
                        feature_meta[penalty_factor_meta_col] == 0
                        or feature_meta['Weight'] != 0]
                else:
                    feature_meta = feature_meta.loc[feature_meta['Weight']
                                                    != 0]
    feature_meta.index.rename('Feature', inplace=True)
    return feature_meta


def add_param_cv_scores(search, param_grid_dict, param_cv_scores=None):
    if param_cv_scores is None:
        param_cv_scores = {}
    for param, param_values in param_grid_dict.items():
        if len(param_values) == 1:
            continue
        param_cv_values = search.cv_results_['param_{}'.format(param)]
        if any(isinstance(v, BaseEstimator) for v in param_cv_values):
            param_cv_values = np.array(
                ['.'.join([type(v).__module__, type(v).__qualname__])
                 if isinstance(v, BaseEstimator) else v
                 for v in param_cv_values])
        if param not in param_cv_scores:
            param_cv_scores[param] = {}
        for metric in metrics:
            if metric not in param_cv_scores[param]:
                param_cv_scores[param][metric] = {'scores': [], 'stdev': []}
            param_metric_scores = param_cv_scores[param][metric]['scores']
            param_metric_stdev = param_cv_scores[param][metric]['stdev']
            for param_value_idx, param_value in enumerate(param_values):
                mean_cv_scores = (search.cv_results_
                                  ['mean_test_{}'.format(metric)]
                                  [param_cv_values == param_value])
                std_cv_scores = (search.cv_results_
                                 ['std_test_{}'.format(metric)]
                                 [param_cv_values == param_value])
                if mean_cv_scores.size > 0:
                    if param_value_idx < len(param_metric_scores):
                        param_metric_scores[param_value_idx] = np.append(
                            param_metric_scores[param_value_idx],
                            mean_cv_scores[np.argmax(mean_cv_scores)])
                        param_metric_stdev[param_value_idx] = np.append(
                            param_metric_stdev[param_value_idx],
                            std_cv_scores[np.argmax(mean_cv_scores)])
                    else:
                        param_metric_scores.append(np.array(
                            [mean_cv_scores[np.argmax(mean_cv_scores)]]))
                        param_metric_stdev.append(np.array(
                            [std_cv_scores[np.argmax(mean_cv_scores)]]))
                elif param_value_idx < len(param_metric_scores):
                    param_metric_scores[param_value_idx] = np.append(
                        param_metric_scores[param_value_idx], [np.nan])
                    param_metric_stdev[param_value_idx] = np.append(
                        param_metric_stdev[param_value_idx], [np.nan])
                else:
                    param_metric_scores.append(np.array([np.nan]))
                    param_metric_stdev.append(np.array([np.nan]))
    return param_cv_scores


def plot_param_cv_metrics(model_name, param_grid_dict, param_cv_scores):
    metric_colors = sns.color_palette(args.sns_color_palette, len(metrics))
    for param in param_cv_scores:
        mean_cv_scores, std_cv_scores = {}, {}
        for metric in metrics:
            param_metric_scores = param_cv_scores[param][metric]['scores']
            param_metric_stdev = param_cv_scores[param][metric]['stdev']
            if any(len(scores) > 1 for scores in param_metric_scores):
                mean_cv_scores[metric], std_cv_scores[metric] = [], []
                for param_value_scores in param_metric_scores:
                    mean_cv_scores[metric].append(
                        np.nanmean(param_value_scores))
                    std_cv_scores[metric].append(
                        np.nanstd(param_value_scores))
            else:
                mean_cv_scores[metric] = np.ravel(param_metric_scores)
                std_cv_scores[metric] = np.ravel(param_metric_stdev)
        plt.figure(figsize=(args.fig_width, args.fig_height))
        param_type = get_param_type(param)
        if param_type in params_lin_xticks:
            x_axis = param_grid_dict[param]
            if all(0 <= x <= 1 for x in x_axis):
                if len(x_axis) <= 15:
                    plt.xticks(x_axis)
            elif len(x_axis) <= 30:
                plt.xticks(x_axis)
        elif param_type in params_log_xticks:
            x_axis = np.ravel(param_grid_dict[param])
            plt.xscale('log', base=(2 if np.all(np.frexp(x_axis)[0] == 0.5)
                                    else 10))
        elif param_type in params_fixed_xticks:
            x_axis = range(len(param_grid_dict[param]))
            xtick_labels = [v.split('.')[-1]
                            if param_type in pipe_step_types
                            and not args.long_label_names
                            and v is not None else str(v)
                            for v in param_grid_dict[param]]
            plt.xticks(x_axis, xtick_labels)
        else:
            raise RuntimeError('No ticks config exists for {}'
                               .format(param_type))
        plt.xlim([min(x_axis), max(x_axis)])
        plt.title('Effect of {} on CV Performance Metrics\n{}'
                  .format(param, model_name), fontsize=args.title_font_size)
        plt.xlabel(param, fontsize=args.axis_font_size)
        plt.ylabel('CV Score', fontsize=args.axis_font_size)
        for metric_idx, metric in enumerate(metrics):
            plt.plot(x_axis, mean_cv_scores[metric],
                     color=metric_colors[metric_idx], lw=2, alpha=0.8,
                     label='Mean {}'.format(metric_label[metric]))
            plt.fill_between(x_axis,
                             [m - s for m, s in zip(mean_cv_scores[metric],
                                                    std_cv_scores[metric])],
                             [m + s for m, s in zip(mean_cv_scores[metric],
                                                    std_cv_scores[metric])],
                             alpha=0.1, color=metric_colors[metric_idx],
                             label=(r'$\pm$ 1 std. dev.'
                                    if metric_idx == len(metrics) - 1
                                    else None))
        plt.legend(loc='lower right', fontsize='medium')
        plt.tick_params(labelsize=args.axis_font_size)
        plt.grid(True, alpha=0.3)


def get_coxnet_max_num_alphas(search):
    param_combos = ParameterGrid(search.param_grid)
    max_num_alphas = 0
    pipe = search.estimator
    srv_step_name = pipe.steps[-1][0]
    cnet_srv_n_param = '{}__estimator__n_alphas'.format(srv_step_name)
    for params in param_combos:
        if (isinstance(pipe[-1], MetaCoxnetSurvivalAnalysis)
                or (srv_step_name in params and isinstance(
                    params[srv_step_name], MetaCoxnetSurvivalAnalysis))):
            max_num_alphas = max(max_num_alphas,
                                 params[cnet_srv_n_param]
                                 if cnet_srv_n_param in params else
                                 params[srv_step_name].estimator.n_alphas
                                 if srv_step_name in params else
                                 pipe[-1].estimator.n_alphas)
    return max_num_alphas


def add_coxnet_alpha_param_grid(search, X, y, pipe_fit_params):
    cnet_pipes = []
    param_combos = ParameterGrid(search.param_grid)
    pipe = search.estimator
    srv_step_name = pipe.steps[-1][0]
    for params in param_combos:
        if (isinstance(pipe[-1], MetaCoxnetSurvivalAnalysis)
                or (srv_step_name in params and isinstance(
                    params[srv_step_name], MetaCoxnetSurvivalAnalysis))):
            cnet_pipe = clone(pipe)
            cnet_pipe.set_params(**params)
            cnet_pipe.steps[-1] = (srv_step_name, cnet_pipe[-1].estimator)
            for param in cnet_pipe.get_params(deep=True).keys():
                param_parts = param.split('__')
                if param_parts[-1] == 'fit_baseline_model':
                    cnet_pipe.set_params(**{param: False})
            cnet_pipes.append(cnet_pipe)
    print('Generating CoxnetSurvivalAnalysis alpha path for {} pipeline{}'
          .format(len(cnet_pipes), 's' if len(cnet_pipes) > 1 else ''),
          flush=True, end='\n' if args.scv_verbose > 0 else ' ')
    fitted_cnet_pipes = Parallel(
        backend=args.parallel_backend, n_jobs=args.n_jobs,
        verbose=args.scv_verbose)(
            delayed(fit_pipeline)(X, y, cnet_pipe.steps, params=None,
                                  param_routing=cnet_pipe.param_routing,
                                  fit_params=pipe_fit_params)
            for cnet_pipe in cnet_pipes)
    if args.scv_verbose == 0:
        print(flush=True)
    if all(p is None for p in fitted_cnet_pipes):
        raise RuntimeError('All CoxnetSurvivalAnalysis alpha path pipelines '
                           'failed')
    param_grid = []
    cnet_pipes_idx = 0
    cnet_srv_a_param = '{}__alpha'.format(srv_step_name)
    for params in param_combos:
        param_grid.append({k: [v] for k, v in params.items()})
        if (isinstance(pipe[-1], MetaCoxnetSurvivalAnalysis)
                or (srv_step_name in params and isinstance(
                    params[srv_step_name], MetaCoxnetSurvivalAnalysis))):
            if fitted_cnet_pipes[cnet_pipes_idx] is not None:
                param_grid[-1][cnet_srv_a_param] = (
                    fitted_cnet_pipes[cnet_pipes_idx][-1].alphas_)
            else:
                del param_grid[-1]
            cnet_pipes_idx += 1
    search.set_params(param_grid=param_grid)
    if args.verbose > 1:
        print('Param grid:')
        pprint(param_grid)
    return search


def update_coxnet_param_grid_dict(search, param_grid_dict):
    pipe = search.estimator
    srv_step_name = pipe.steps[-1][0]
    cnet_srv_a_param = '{}__alpha'.format(srv_step_name)
    cnet_srv_l_param = '{}__estimator__l1_ratio'.format(srv_step_name)
    cnet_srv_n_param = '{}__estimator__n_alphas'.format(srv_step_name)
    if any(p in search.best_params_ for p in (cnet_srv_l_param,
                                              cnet_srv_n_param)):
        best_alpha_condition = {k: v for k, v in search.best_params_.items()
                                if k in (cnet_srv_l_param, cnet_srv_n_param)}
        param_grid_dict[cnet_srv_a_param] = list(filter(
            lambda params: all(params[k] == [v] for k, v in
                               best_alpha_condition.items()),
            search.param_grid))[0][cnet_srv_a_param]
    else:
        param_grid_dict[cnet_srv_a_param] = (
            search.param_grid[0][cnet_srv_a_param])
    return param_grid_dict


def unset_pipe_memory(pipe):
    for param, param_value in pipe.get_params(deep=True).items():
        if isinstance(param_value, Memory):
            pipe.set_params(**{param: None})
    if (isinstance(pipe[0], ColumnTransformer)
            and hasattr(pipe[0], 'transformers_')):
        for _, trf_transformer, _ in pipe[0].transformers_:
            if isinstance(trf_transformer, Pipeline):
                unset_pipe_memory(trf_transformer)
    return pipe


def run_model():
    (dataset_name, X, y, groups, group_weights, sample_weights, sample_meta,
     feature_meta) = load_dataset(args.dataset)
    pipe, param_grid, param_grid_dict = setup_pipe_and_param_grid(X)
    pipe_has_penalty_factor = False
    for param in pipe.get_params(deep=True).keys():
        param_parts = param.split('__')
        if param_parts[-1] == 'penalty_factor_meta_col':
            pipe.set_params(**{param: penalty_factor_meta_col})
            pipe_has_penalty_factor = True
    for params in param_grid:
        for param_values in params.values():
            if any(isinstance(v, BaseEstimator) for v in param_values):
                for estimator in param_values:
                    for param in estimator.get_params(deep=True).keys():
                        param_parts = param.split('__')
                        if param_parts[-1] == 'penalty_factor_meta_col':
                            estimator.set_params(
                                **{param: penalty_factor_meta_col})
                            pipe_has_penalty_factor = True
    if not pipe_has_penalty_factor:
        feature_meta.drop(columns=[penalty_factor_meta_col],
                          inplace=True)
    if groups is not None:
        search_param_routing = {'cv': 'groups', 'estimator': [], 'scoring': []}
    else:
        search_param_routing = None
    if pipe.param_routing:
        if search_param_routing is None:
            search_param_routing = {'estimator': [], 'scoring': []}
        for param in [p for l in pipe.param_routing.values() for p in l]:
            if param not in search_param_routing['estimator']:
                search_param_routing['estimator'].append(param)
                search_param_routing['scoring'].append(param)
    test_split_params = {'groups': groups} if groups is not None else {}
    if analysis == 'surv':
        if groups is None:
            test_splitter = SurvivalStratifiedShuffleSplit(
                n_splits=test_splits, test_size=test_size,
                random_state=random_seed)
        else:
            test_splitter = SurvivalStratifiedSampleFromGroupShuffleSplit(
                n_splits=test_splits, test_size=test_size,
                random_state=random_seed)
            test_split_params['weights'] = group_weights
        cv_splitter = RepeatedSurvivalStratifiedKFold(
            n_splits=scv_splits, n_repeats=scv_repeats,
            random_state=random_seed)
    elif (groups is None or 'sample_weight'
          not in search_param_routing['estimator']):
        test_splitter = RepeatedStratifiedKFold(
            n_splits=test_splits, n_repeats=test_repeats,
            random_state=random_seed)
        cv_splitter = RepeatedStratifiedKFold(
            n_splits=scv_splits, n_repeats=scv_repeats,
            random_state=random_seed)
    else:
        test_splitter = RepeatedStratifiedGroupKFold(
            n_splits=test_splits, n_repeats=test_repeats,
            random_state=random_seed)
        cv_splitter = RepeatedStratifiedGroupKFold(
            n_splits=scv_splits, n_repeats=scv_repeats,
            random_state=random_seed)
    if refit_metric == 'score':
        scv_scoring = None
        scv_refit = True
    else:
        scv_scoring = metrics
        scv_refit = refit_metric
    scv_error_score = 0 if args.model_type in ('cnet', 'lgr') else 'raise'
    base_search = ExtendedGridSearchCV(
        pipe, cv=cv_splitter, error_score=scv_error_score,
        n_jobs=None, param_grid=param_grid, param_routing=search_param_routing,
        refit=scv_refit, return_train_score=False, scoring=scv_scoring,
        verbose=args.scv_verbose)
    if args.verbose > 0:
        print(base_search.__repr__(N_CHAR_MAX=10000))
        if param_grid_dict:
            print('Param grid dict:')
            pprint(param_grid_dict)
    if args.verbose > 0 or args.scv_verbose > 0:
        print('Dataset:', dataset_name, X.shape, end=' ')
        if isinstance(pipe[0], ColumnTransformer):
            print(col_trf_info(pipe[0]))
        else:
            print()
    if args.verbose > 0:
        if groups is not None:
            print('Groups:')
            pprint(groups)
            if group_weights is not None:
                print('Group weights:')
                pprint(group_weights)
        if (sample_weights is not None and 'sample_weight' in
                search_param_routing['estimator']):
            print('Sample weights:')
            pprint(sample_weights)
        print('Test CV:', end=' ')
        pprint(test_splitter)
    model_name = '_'.join([dataset_name.rpartition('_')[0], args.model_type])
    if args.load_only:
        sys.exit()
    pipe_fit_params = {}
    if search_param_routing:
        if 'sample_meta' in search_param_routing['estimator']:
            pipe_fit_params['sample_meta'] = sample_meta
        if 'feature_meta' in search_param_routing['estimator']:
            pipe_fit_params['feature_meta'] = feature_meta
        if 'sample_weight' in search_param_routing['estimator']:
            pipe_fit_params['sample_weight'] = sample_weights
    search_fit_params = pipe_fit_params.copy()
    if groups is not None:
        search_fit_params['groups'] = groups
    split_models = []
    split_results = []
    param_cv_scores = {}
    if analysis == 'surv':
        param_grid_dict_alphas = None
        risk_scores = None
        srv_step_name = pipe.steps[-1][0]
        cnet_srv_a_param = '{}__alpha'.format(srv_step_name)
        cnet_srv_max_num_alphas = get_coxnet_max_num_alphas(base_search)
    else:
        print('Generating permutation test input data')
        random_state = check_random_state(random_seed)
        perm_ys, perm_split_idxs = zip(*Parallel(
            backend=args.parallel_backend, n_jobs=args.n_jobs,
            verbose=args.perm_verbose)(delayed(get_perm_test_split_data)(
                X, shuffle_y(y, groups, random_state), test_splitter,
                cv_params=test_split_params) for _ in range(args.n_perms)))
        split_perm_idxs = [*zip(*perm_split_idxs)]
    for split_idx, (train_idxs, test_idxs) in enumerate(
            test_splitter.split(X, y, **test_split_params)):
        split_pipe_fit_params = {
            k: (v.iloc[train_idxs] if k in ('sample_meta')
                else v[train_idxs] if k in ('sample_weight')
                else v)
            for k, v in pipe_fit_params.items() if v is not None}
        split_search_fit_params = split_pipe_fit_params.copy()
        if groups is not None:
            split_search_fit_params['groups'] = groups[train_idxs]
        if analysis == 'surv':
            search = add_coxnet_alpha_param_grid(
                clone(base_search), X.iloc[train_idxs], y[train_idxs],
                split_pipe_fit_params)
        else:
            search = clone(base_search)
        try:
            with parallel_backend(args.parallel_backend, n_jobs=args.n_jobs,
                                  inner_max_num_threads=1):
                search.fit(X.iloc[train_idxs], y[train_idxs],
                           **split_search_fit_params)
            best_index = search.best_index_
            best_params = search.best_params_
            best_pipe = search.best_estimator_
            split_scores = {'cv': {}}
            for metric in metrics:
                split_scores['cv'][metric] = search.cv_results_[
                    'mean_test_{}'.format(metric)][best_index]
            split_pipe_predict_params = {
                k: v.iloc[test_idxs] if k in ('sample_meta') else v
                for k, v in pipe_fit_params.items()
                if k not in ('sample_weight') and v is not None}
            split_score_params = {
                'sample_weight': (sample_weights[test_idxs]
                                  if sample_weights is not None else None)}
            split_scores['te'] = calculate_test_scores(
                best_pipe, X.iloc[test_idxs], y[test_idxs],
                metrics, predict_params=split_pipe_predict_params,
                score_params=split_score_params)
            if analysis == 'surv':
                surv_funcs = best_pipe.predict_survival_function(
                    X.iloc[test_idxs], **split_pipe_predict_params)
            else:
                print('Running permutation test ({:d} permutations)'
                      .format(args.n_perms))
                split_perm_scores = Parallel(
                    backend=args.parallel_backend, n_jobs=args.n_jobs,
                    verbose=args.perm_verbose)(
                        delayed(fit_and_score)(
                            unset_pipe_memory(clone(best_pipe)),
                            X.iloc[perm_train_idxs],
                            perm_y[perm_train_idxs],
                            X.iloc[perm_test_idxs],
                            perm_y[perm_test_idxs],
                            refit_metric,
                            fit_params={
                                k: (v.iloc[perm_train_idxs]
                                    if k in ('sample_meta')
                                    else v[perm_train_idxs]
                                    if k in ('sample_weight')
                                    else v)
                                for k, v in pipe_fit_params.items()
                                if v is not None},
                            predict_params={
                                k: (v.iloc[perm_test_idxs]
                                    if k in ('sample_meta')
                                    else v)
                                for k, v in pipe_fit_params.items()
                                if k not in ('sample_weight')
                                and v is not None},
                            score_params={
                                'sample_weight': (
                                    sample_weights[perm_test_idxs]
                                    if sample_weights is not None
                                    else None)})
                        for perm_y, (perm_train_idxs, perm_test_idxs) in (
                            zip(perm_ys, split_perm_idxs[split_idx])))
        except Exception as e:
            if search.error_score == 'raise':
                raise
            if args.verbose > 0:
                print('Dataset:', dataset_name, ' Split: {:>{width}d}'
                      .format(split_idx + 1, width=len(str(test_splits))),
                      end=' ', flush=True)
            warnings.formatwarning = warning_format
            warnings.warn('Estimator fit/scoring failed. This outer CV '
                          'train-test split will be ignored. Details: {}'
                          .format(format_exception_only(type(e), e)[0]),
                          category=FitFailedWarning)
            best_pipe = None
            split_result = None
        else:
            if analysis == 'surv':
                param_grid_dict = update_coxnet_param_grid_dict(
                    search, param_grid_dict)
                full_alpha_path = np.logspace(
                    np.log10(param_grid_dict[cnet_srv_a_param][0]),
                    np.log10(param_grid_dict[cnet_srv_a_param][-1]),
                    num=cnet_srv_max_num_alphas)
                param_grid_dict_alphas = (
                    np.vstack((param_grid_dict_alphas, full_alpha_path))
                    if param_grid_dict_alphas is not None else
                    full_alpha_path)
            param_cv_scores = add_param_cv_scores(search, param_grid_dict,
                                                  param_cv_scores)
            final_feature_meta = get_final_feature_meta(best_pipe,
                                                        feature_meta)
            if args.verbose > 0:
                print('Model:', model_name, ' Split: {:>{width}d}'
                      .format(split_idx + 1, width=len(str(test_splits))),
                      end=' ')
                for metric in metrics:
                    print(' {} (CV / Test): {:.4f} / {:.4f}'.format(
                        metric_label[metric], split_scores['cv'][metric],
                        split_scores['te'][metric]), end=' ')
                    if metric == 'average_precision':
                        print(' PR AUC Test: {:.4f}'.format(
                            split_scores['te']['pr_auc']), end=' ')
                print(' Params:', {
                    k: ('.'.join([type(v).__module__,
                                  type(v).__qualname__])
                        if isinstance(v, BaseEstimator) else v)
                    for k, v in best_params.items()}, end=' ')
            if penalty_factor_meta_col in final_feature_meta.columns:
                num_features = final_feature_meta.loc[
                    final_feature_meta[penalty_factor_meta_col] != 0].shape[0]
            else:
                num_features = final_feature_meta.shape[0]
            print(' Features: {:.0f}'.format(num_features))
            if args.verbose > 1:
                if 'Weight' in final_feature_meta.columns:
                    print(tabulate(final_feature_meta.iloc[
                        (-final_feature_meta['Weight'].abs()).argsort()],
                                   floatfmt='.6e', headers='keys'))
                else:
                    print(tabulate(final_feature_meta, headers='keys'))
            split_result = {'feature_meta': final_feature_meta,
                            'scores': split_scores}
            if analysis == 'surv':
                split_result['surv_funcs'] = surv_funcs
                split_risk_score_col = 'Risk Score {:d}'.format(split_idx + 1)
                split_risk_scores = pd.DataFrame(
                    {split_risk_score_col: split_scores['te']['y_pred']},
                    index=sample_meta.index[test_idxs])
                risk_scores = (risk_scores.join(split_risk_scores, how='outer')
                               if risk_scores is not None else
                               split_risk_scores)
            else:
                split_result['perm_scores'] = split_perm_scores
        split_results.append(split_result)
        if best_pipe is not None:
            best_pipe = unset_pipe_memory(best_pipe)
        split_models.append(best_pipe)
        memory.clear(warn=False)
    results_dir = '{}/{}'.format(out_dir, model_name)
    os.makedirs(results_dir, mode=0o755, exist_ok=True)
    dump(split_models, '{}/{}_split_models.pkl'
         .format(results_dir, model_name))
    dump(split_results, '{}/{}_split_results.pkl'
         .format(results_dir, model_name))
    dump(param_cv_scores, '{}/{}_param_cv_scores.pkl'
         .format(results_dir, model_name))
    if analysis == 'surv':
        dump(risk_scores, '{}/{}_risk_scores.pkl'
             .format(results_dir, model_name))
        r_base.saveRDS(risk_scores, '{}/{}_risk_scores.rds'
                       .format(results_dir, model_name))
    if analysis == 'surv':
        param_grid_dict[cnet_srv_a_param] = np.mean(
            param_grid_dict_alphas, axis=0)
    scores = {'cv': {}, 'te': {}}
    num_features = []
    if analysis == 'resp':
        perm_scores = []
    for split_result in split_results:
        if split_result is None:
            continue
        for metric in metrics:
            if metric not in scores['cv']:
                scores['cv'][metric] = []
                scores['te'][metric] = []
            scores['cv'][metric].append(
                split_result['scores']['cv'][metric])
            scores['te'][metric].append(
                split_result['scores']['te'][metric])
            if metric == 'average_precision':
                if 'pr_auc' not in scores['te']:
                    scores['te']['pr_auc'] = []
                scores['te']['pr_auc'].append(
                    split_result['scores']['te']['pr_auc'])
        split_feature_meta = split_result['feature_meta']
        if penalty_factor_meta_col in split_feature_meta.columns:
            num_features.append(split_feature_meta.loc[
                split_feature_meta[penalty_factor_meta_col]
                != 0].shape[0])
        else:
            num_features.append(split_feature_meta.shape[0])
        if analysis == 'resp':
            perm_scores.append(split_result['perm_scores'])
    if analysis == 'resp':
        perm_scores = np.mean(perm_scores, axis=0)
        true_score = np.mean(scores['te'][refit_metric])
        perm_pvalue = ((np.sum(perm_scores >= true_score) + 1.0)
                       / (args.n_perms + 1))
        perm_results = {'true_score': true_score,
                        'scores': perm_scores,
                        'pvalue': perm_pvalue}
        dump(perm_results, '{}/{}_perm_results.pkl'
             .format(results_dir, model_name))
    print('Model:', model_name, end=' ')
    for metric in metrics:
        print(' Mean {} (CV / Test): {:.4f} / {:.4f}'.format(
            metric_label[metric], np.mean(scores['cv'][metric]),
            np.mean(scores['te'][metric])), end=' ')
        if metric == 'average_precision':
            print(' Mean PR AUC Test: {:.4f}'.format(
                np.mean(scores['te']['pr_auc'])), end=' ')
    print(' Mean Features: {:.0f}'.format(np.mean(num_features)), end=' ')
    if analysis == 'resp':
        print(' Permutation Test: True {} = {:.4f} p = {:.4f}'.format(
            metric_label[refit_metric], true_score, perm_pvalue))
    else:
        print()
    # feature mean rankings and scores
    feature_annots = None
    feature_weights = None
    feature_scores = {}
    for split_idx, split_result in enumerate(split_results):
        if split_result is None:
            continue
        split_feature_meta = split_result['feature_meta']
        if feature_meta.columns.any():
            if feature_annots is None:
                feature_annots = split_feature_meta[feature_meta.columns]
            else:
                feature_annots = pd.concat(
                    [feature_annots,
                     split_feature_meta[feature_meta.columns]], axis=0)
        elif feature_annots is None:
            feature_annots = pd.DataFrame(index=split_feature_meta.index)
        else:
            feature_annots = pd.concat(
                [feature_annots,
                 pd.DataFrame(index=split_feature_meta.index)], axis=0)
        if 'Weight' in split_feature_meta.columns:
            if feature_weights is None:
                feature_weights = split_feature_meta[['Weight']].copy()
            else:
                feature_weights = feature_weights.join(
                    split_feature_meta[['Weight']], how='outer')
            feature_weights.rename(
                columns={'Weight': 'Weight {:d}'.format(split_idx + 1)},
                inplace=True)
        for metric in metrics:
            if metric not in feature_scores:
                feature_scores[metric] = pd.DataFrame(
                    split_result['scores']['te'][metric], columns=[metric],
                    index=split_feature_meta.index)
            else:
                feature_scores[metric] = feature_scores[metric].join(
                    pd.DataFrame(split_result['scores']['te'][metric],
                                 columns=[metric],
                                 index=split_feature_meta.index),
                    how='outer')
            feature_scores[metric].rename(columns={metric: split_idx},
                                          inplace=True)
    feature_annots = feature_annots.loc[
        ~feature_annots.index.duplicated(keep='first')]
    feature_frequency = None
    feature_results = None
    feature_results_floatfmt = ['']
    if feature_weights is not None:
        feature_ranks = feature_weights.abs().rank(
            ascending=False, method='min', na_option='keep')
        feature_ranks.fillna(feature_ranks.count(axis=0) + 1, inplace=True)
        feature_frequency = feature_weights.count(axis=1)
        feature_weights.fillna(0, inplace=True)
        feature_results = feature_annots.reindex(index=feature_ranks.index,
                                                 fill_value='')
        for feature_annot_col in feature_annots.columns:
            if is_integer_dtype(feature_annots[feature_annot_col]):
                feature_results_floatfmt.append('.0f')
            elif is_float_dtype(feature_annots[feature_annot_col]):
                feature_results_floatfmt.append('.{:d}f'.format(
                    max(abs(Decimal(f).as_tuple().exponent)
                        for f in (feature_annots[feature_annot_col]
                                  .astype(str)))))
            else:
                feature_results_floatfmt.append('')
        feature_results['Frequency'] = feature_frequency
        feature_results['Mean Weight Rank'] = feature_ranks.mean(axis=1)
        feature_results['Mean Weight'] = feature_weights.mean(axis=1)
        feature_results_floatfmt.extend(['.0f', '.1f', '.6e'])
    for metric in metrics:
        if feature_results is None:
            feature_results = feature_annots.reindex(
                index=feature_scores[metric].index, fill_value='')
            for feature_annot_col in feature_annots.columns:
                if is_integer_dtype(feature_annots[feature_annot_col]):
                    feature_results_floatfmt.append('.0f')
                elif is_float_dtype(feature_annots[feature_annot_col]):
                    feature_results_floatfmt.append('.{:d}f'.format(
                        max(abs(Decimal(f).as_tuple().exponent)
                            for f in (feature_annots[feature_annot_col]
                                      .astype(str)))))
                else:
                    feature_results_floatfmt.append('')
            feature_frequency = feature_scores[metric].count(axis=1)
            feature_results['Frequency'] = feature_frequency
            feature_results_floatfmt.append('.0f')
        feature_scores[metric].fillna(0.5, inplace=True)
        if feature_scores[metric].mean(axis=1).nunique() > 1:
            feature_results = feature_results.join(
                pd.DataFrame({
                    'Mean Test {}'.format(metric_label[metric]):
                        feature_scores[metric].mean(axis=1)}),
                how='left')
            feature_results_floatfmt.append('.4f')
    dump(feature_results, '{}/{}_feature_results.pkl'
         .format(results_dir, model_name))
    r_base.saveRDS(feature_results, '{}/{}_feature_results.rds'
                   .format(results_dir, model_name))
    if feature_weights is not None:
        dump(feature_weights, '{}/{}_feature_weights.pkl'
             .format(results_dir, model_name))
        r_base.saveRDS(feature_weights, '{}/{}_feature_weights.rds'
                       .format(results_dir, model_name))
    if args.verbose > 0:
        print('Overall Feature Ranking:')
        if feature_weights is not None:
            print(tabulate(
                feature_results.sort_values(by='Mean Weight Rank'),
                floatfmt=feature_results_floatfmt, headers='keys'))
        else:
            print(tabulate(
                feature_results.sort_values(by='Mean Test {}'.format(
                    metric_label[refit_metric]), ascending=False),
                floatfmt=feature_results_floatfmt, headers='keys'))
    plot_param_cv_metrics(model_name, param_grid_dict, param_cv_scores)
    # plot roc and pr curves
    metric_colors = sns.color_palette(args.sns_color_palette, len(metrics))
    if 'roc_auc' in metrics:
        metric_idx = metrics.index('roc_auc')
        plt.figure(figsize=(args.fig_width, args.fig_height))
        plt.title('ROC Curve\n{}'.format(model_name),
                  fontsize=args.title_font_size)
        plt.xlabel('False Positive Rate', fontsize=args.axis_font_size)
        plt.ylabel('True Positive Rate', fontsize=args.axis_font_size)
        plt.xlim([-0.01, 1.01])
        plt.ylim([-0.01, 1.01])
        tprs = []
        mean_fpr = np.linspace(0, 1, 100)
        for split_result in split_results:
            if split_result is None:
                continue
            tprs.append(np.interp(mean_fpr,
                                  split_result['scores']['te']['fpr'],
                                  split_result['scores']['te']['tpr']))
            tprs[-1][0] = 0.0
            plt.plot(split_result['scores']['te']['fpr'],
                     split_result['scores']['te']['tpr'], alpha=0.2,
                     color='darkgrey', lw=1)
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_roc_auc = np.mean(scores['te']['roc_auc'])
        std_roc_auc = np.std(scores['te']['roc_auc'])
        mean_num_features = np.mean(num_features)
        std_num_features = np.std(num_features)
        plt.plot(mean_fpr, mean_tpr, lw=2, alpha=0.8,
                 color=metric_colors[metric_idx],
                 label=(r'Test Mean ROC (AUC = {:.4f} $\pm$ {:.2f}, '
                        r'Features = {:.0f} $\pm$ {:.0f})').format(
                            mean_roc_auc, std_roc_auc, mean_num_features,
                            std_num_features))
        std_tpr = np.std(tprs, axis=0)
        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        plt.fill_between(mean_fpr, tprs_lower, tprs_upper, alpha=0.1,
                         color='grey', label=r'$\pm$ 1 std. dev.')
        plt.plot([0, 1], [0, 1], alpha=0.2, color='grey',
                 linestyle='--', lw=2, label='Chance')
        plt.legend(loc='lower right', fontsize='medium')
        plt.tick_params(labelsize=args.axis_font_size)
        plt.grid(False)
    if 'average_precision' in metrics:
        metric_idx = metrics.index('average_precision')
        plt.figure(figsize=(args.fig_width, args.fig_height))
        plt.title('PR Curve\n{}'.format(model_name),
                  fontsize=args.title_font_size)
        plt.xlabel('Recall', fontsize=args.axis_font_size)
        plt.ylabel('Precision', fontsize=args.axis_font_size)
        plt.xlim([-0.01, 1.01])
        plt.ylim([-0.01, 1.01])
        pres = []
        mean_rec = np.linspace(0, 1, 100)
        for split_result in split_results:
            if split_result is None:
                continue
            pres.append(np.interp(
                mean_rec, split_result['scores']['te']['rec'][::-1],
                split_result['scores']['te']['pre'][::-1]))
            plt.step(split_result['scores']['te']['rec'],
                     split_result['scores']['te']['pre'], alpha=0.2,
                     color='darkgrey', lw=1, where='post')
        mean_pre = np.mean(pres, axis=0)
        mean_pr_auc = np.mean(scores['te']['pr_auc'])
        std_pr_auc = np.std(scores['te']['pr_auc'])
        mean_num_features = np.mean(num_features)
        std_num_features = np.std(num_features)
        plt.step(mean_rec, mean_pre, alpha=0.8, lw=2,
                 color=metric_colors[metric_idx], where='post',
                 label=(r'Test Mean PR (AUC = {:.4f} $\pm$ {:.2f}, '
                        r'Features = {:.0f} $\pm$ {:.0f})').format(
                            mean_pr_auc, std_pr_auc, mean_num_features,
                            std_num_features))
        std_pre = np.std(pres, axis=0)
        pres_upper = np.minimum(mean_pre + std_pre, 1)
        pres_lower = np.maximum(mean_pre - std_pre, 0)
        plt.fill_between(mean_rec, pres_lower, pres_upper, alpha=0.1,
                         color='grey', label=r'$\pm$ 1 std. dev.')
        plt.legend(loc='lower right', fontsize='medium')
        plt.tick_params(labelsize=args.axis_font_size)
        plt.grid(False)
    # plot permutation test histogram
    if analysis == 'resp':
        metric_idx = metrics.index(refit_metric)
        _, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
        plt.title('Permutation Test\n{}'.format(model_name),
                  fontsize=args.title_font_size)
        # freedman-draconis rule
        bins = round((np.max(perm_scores) - np.min(perm_scores))
                     / (2 * iqr(perm_scores) / np.cbrt(perm_scores.size)))
        sns.histplot(perm_scores, bins=bins, kde=True,
                     color=metric_colors[metric_idx], stat='probability',
                     edgecolor='white')
        plt.axvline(true_score, ls='--', color='darkgrey')
        ax.add_artist(AnchoredText(
            r'True {} = {:.2f}' '\n' r'$\itp$ = {:.2e}'
            .format(metric_label[refit_metric], true_score, perm_pvalue),
            loc='upper left', frameon=False,
            prop={'size': args.axis_font_size}))
        plt.xticks(np.arange(0.0, 1.1, 0.2))
        plt.xlabel(metric_label[refit_metric], fontsize=args.axis_font_size)
        plt.ylabel('Probability', fontsize=args.axis_font_size)
        plt.tick_params(labelsize=args.axis_font_size)
    for fig_num in plt.get_fignums():
        plt.figure(fig_num, constrained_layout=True)
        for fig_fmt in args.fig_format:
            plt.savefig('{}/Figure_{:d}.{}'.format(results_dir, fig_num,
                                                   fig_fmt),
                        bbox_inches='tight', format=fig_fmt)


def run_cleanup():
    rmtree(cachedir)
    for rtmp_dir in glob('{}/Rtmp*/'.format(args.tmp_dir)):
        rmtree(rtmp_dir)


def dir_path(path):
    os.makedirs(path, mode=0o755, exist_ok=True)
    return path


parser = ArgumentParser()
parser.add_argument('--dataset', type=str, required=True, help='dataset')
parser.add_argument('--model-type', type=str, required=True,
                    choices=['cnet', 'rfe', 'lgr', 'edger', 'limma'],
                    help='model type')
parser.add_argument('--scv-splits', type=int, help='num inner cv splits')
parser.add_argument('--scv-repeats', type=int, help='num inner cv repeats')
parser.add_argument('--test-splits', type=int, help='num outer test splits')
parser.add_argument('--test-repeats', type=int, help='num outer test repeats')
parser.add_argument('--test-size', type=float, help='outer cv test size')
parser.add_argument('--scv-verbose', type=int, default=0, help='scv verbosity')
parser.add_argument('--n-jobs', type=int, default=-2, help='num parallel jobs')
parser.add_argument('--parallel-backend', type=str, default='loky',
                    help='joblib parallel backend')
parser.add_argument('--n-perms', type=int, default=1000,
                    help='permutation test n permutations')
parser.add_argument('--perm-verbose', type=int, default=0,
                    help='permutation test verbosity')
parser.add_argument('--sns-color-palette', type=str, default='hls',
                    help='Seaborn/matplotlib color palette')
parser.add_argument('--title-font-size', type=int, default=14,
                    help='figure title font size')
parser.add_argument('--axis-font-size', type=int, default=14,
                    help='figure axis font size')
parser.add_argument('--long-label-names', default=False, action='store_true',
                    help='figure long label names')
parser.add_argument('--fig-width', type=float, default=10,
                    help='figure width')
parser.add_argument('--fig-height', type=float, default=10,
                    help='figure height')
parser.add_argument('--fig-format', type=str, nargs='+',
                    choices=['png', 'pdf', 'svg', 'tif'], default=['png'],
                    help='figure format')
parser.add_argument('--results-dir', type=str, default='results/models',
                    help='results dir')
parser.add_argument('--tmp-dir', type=dir_path, default=gettempdir(),
                    help='tmp dir')
parser.add_argument('--verbose', type=int, default=1, help='program verbosity')
parser.add_argument('--load-only', default=False, action='store_true',
                    help='set up model selection and load dataset only')
args = parser.parse_args()

file_basename = os.path.splitext(os.path.split(args.dataset)[1])[0]
_, cancer, analysis, target, data_type, *rest = file_basename.split('_')
if args.model_type in ('edger', 'limma'):
    args.model_type = 'edger' if data_type == 'htseq' else 'limma'
if analysis == 'surv' and args.model_type != 'cnet':
    args.model_type = 'cnet'

out_dir = '{}/{}'.format(args.results_dir, analysis)
os.makedirs(out_dir, mode=0o755, exist_ok=True)

cancer_target = '_'.join([cancer, target])
if analysis == 'surv':
    metrics = ['score']
    test_splits = 100 if args.test_splits is None else args.test_splits
    test_size = 0.25 if args.test_size is None else args.test_size
    scv_repeats = 5 if args.scv_repeats is None else args.scv_repeats
    if args.scv_splits is None:
        scv_splits = (
            2 if cancer_target in ('dlbc_os', 'pcpg_os', 'tgct_os') else
            3 if cancer_target in ('chol_os', 'chol_pfi', 'dlbc_pfi',
                                   'kich_os', 'thym_os') else 4)
    else:
        scv_splits = args.scv_splits
else:
    metrics = ['roc_auc', 'balanced_accuracy', 'average_precision']
    scv_splits = 3 if args.scv_splits is None else args.scv_splits
    scv_repeats = 5 if args.scv_repeats is None else args.scv_repeats
    if args.test_splits is None:
        test_splits = 3 if cancer_target == 'stad_oxaliplatin' else 4
    else:
        test_splits = args.test_splits
    if args.test_repeats is None:
        test_repeats = 33 if cancer_target == 'stad_oxaliplatin' else 25
    else:
        test_repeats = args.test_repeats

random_seed = 777
refit_metric = metrics[0]

r_base = importr('base')
r_biobase = importr('Biobase')
robjects.r('set.seed({:d})'.format(random_seed))

atexit.register(run_cleanup)

cachedir = mkdtemp(dir=args.tmp_dir)
memory = Memory(location=cachedir, verbose=0)

pipe_step_types = ('slr', 'trf', 'clf', 'srv')
params_lin_xticks = [
    'clf__n_features_to_select',
    'slr__estimator__l1_ratio',
    'slr__k',
    'srv__estimator__l1_ratio']
params_log_xticks = [
    'clf__C',
    'clf__estimator__C',
    'slr__estimator__C',
    'srv__alpha']
params_fixed_xticks = []
metric_label = {
    'score': 'C-index',
    'roc_auc': 'ROC AUC',
    'balanced_accuracy': 'BCR',
    'average_precision': 'AVG PRE'}
penalty_factor_meta_col = 'Penalty Factor'
sample_meta_cols = ['gender', 'age_at_diagnosis', 'tumor_stage']
ordinal_encoder_categories = {
    'tumor_stage': ['0', 'i', 'i or ii', 'ii', 'NA', 'iii', 'iv']}

run_model()
