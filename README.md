# tcga-microbiome-prediction

Hermida LC, Gertz EM, and Ruppin E. Analyzing the tumor microbiome to predict
cancer patient survival and drug response.

### Installation

Install and set up [Miniconda3](https://docs.conda.io/en/latest/miniconda.html)

Set up repository and conda environment:

```bash
git clone git@github.com:ruppinlab/tcga-microbiome-prediction.git
cd tcga-microbiome-prediction
conda env create -f envs/tcga-microbiome-prediction.yml
conda activate tcga-microbiome-prediction
```

### Data preprocessing

Download NCI GDC GENCODE v22 GTF file and update Ensembl gene symbols:

```bash
Rscript get_gtf_latest_ensg_annots.R
```

Download and process Poore et al. data and NCI GDC case metadata:

```bash
Rscript process_knight_data_gdc_meta.R
```

Process TCGA survival and drug response phenotypic data:

```bash
Rscript process_surv_resp_pdata.R
```

Download NCI GDC gene expression data and create microbial, gene expression,
and combination ExpressionSet objects (takes ~3GB space)

```bash
Rscript create_esets.R
```

### Clinical covariate Cox and linear SVM models

Run the clinical covariate models and save the results:

```bash
python run_cox_clinical_models.py
python run_svm_clinical_models.py
```

### Coxnet and SVM-RFE models

This should be run on a cluster as there are 282 total models and this is a
very compute-intensive procedure. See the provided Slurm scripts for more
information on porting to other cluster software systems.

To submit and run all the models on a Slurm cluster and save the results:

```bash
python submit_slurm_models.py
```

You can also run individual models and save the results.

Survival:

```bash
./run_model.sh --dataset data/tcga_acc_surv_os_kraken_eset.rds
```

Drug response:

```bash
./run_model.sh --dataset data/tcga_blca_resp_cisplatin_kraken_eset.rds
```

### Model results

Dump model scores to pandas and R dataframes:

```bash
python dump_model_scores.py --model-code cnet
python dump_model_scores.py --model-code rfe
```

You can generate a summary of the model results as a tsv file:

```bash
python summarize_model_results.py
```

### Figures

Time-dependent cumulative/dynamic AUC plots:

```bash
python generate_surv_td_auc_plots.py
```

ROC and PR curve plots:

```bash
python generate_resp_roc_pr_plots.py
```

### Analysis

If all R libraries are installed, and results are in `results`, and
particularly model scores are in `results/model_scores`, with filenames

    cnet_model_scores.rds
    cox_covariate_model_scores.rds
    rfe_model_scores.rds
    svm_covariate_model_scores.rds

then type `make -f Makefile.figures`.

The figures will be generated in the subdirectory `figures`.
