# tcga-microbiome-prediction

Hermida LC, Gertz EM, and Ruppin E. Predicting cancer prognosis and drug
response from the tumor microbiome.

To reproduce the work associated with this project, please follow the steps
below in order.

### Installation

Install and set up [Miniconda3](https://docs.conda.io/en/latest/miniconda.html)

Miniconda3 is designed to be small, and the installation will take 1-5
minutes on a typical computer, depending on the internet connection.

Set up repository and conda environment:

```bash
git clone https://github.com/ruppinlab/tcga-microbiome-prediction.git
cd tcga-microbiome-prediction
conda env create -f envs/tcga-microbiome-prediction.yml
conda activate tcga-microbiome-prediction
```

This time to complete this step is dependent on the internet
connection, but with a typical computer and internet connection takes 1-2
minutes.

### Data preprocessing

Download NCI GDC GENCODE v22 GTF file and update Ensembl gene symbols:

```bash
Rscript get_gtf_ensg_annots.R
```

Download and process Poore et al. microbial abundance data and NCI GDC case
metadata:

```bash
Rscript process_knight_data_gdc_meta.R
```

Process TCGA prognosis and drug response phenotypic data:

```bash
Rscript process_surv_resp_pdata.R
```

Download NCI GDC gene expression data and create microbial abundance, gene
expression, and combination ExpressionSet objects (takes ~3GB space)

```bash
Rscript create_esets.R
```

The most significant time in the data preprocessing step is typically
in the `create_esets.R` script. The actual time required will
depend not only on the speed of the internet connection, but on how
busy the NCI GDC servers are. With a typical internet connection speeds it
takes ~1 hour.

### Clinical covariate models

Run the clinical covariate-only models and save the results:

```bash
python run_surv_clinical_models.py
python run_resp_clinical_models.py
```

### Microbial abundance, gene expression, and combination data type models

Building the ML models for this project was done on the NIH Biowulf
cluster. We recommend running on a cluster of computers because there
are 462 total models and this is a very compute-intensive procedure.
The Biowulf cluster uses the Slurm queuing system software, which is a
typical scheduling system for compute clusters. We've provided our
Slurm scripts, which should be straightforward to port to another job
queuing system.

To submit and run all the models on a Slurm cluster and save the results:

```bash
python submit_slurm_models.py
```

You can also run individual models locally and save the results using
`run_model.sh`. Please note that even running a single model is a
compute-intensive procedure because each "model" is actually comprised of 100
model instances/iterations generated from randomly shuffled train/test data
splits, and each of the 100 model instances undergoes an exhaustive grid search
to tune model hyperparameters and select the best model using cross validation
of randomly shuffled train/validation data splits. So, depending on the data
type, model type, and your available CPU resources, each "model" can take from
a couple hours to days. For example:

Prognosis:

```bash
./run_model.sh --model-type cnet --dataset data/tcga_acc_surv_os_kraken_eset.rds
```

Drug response:

```bash
./run_model.sh --model-type rfe --dataset data/tcga_blca_resp_cisplatin_kraken_eset.rds
```

### Model results

Dump model scores to pandas and R dataframes:

```bash
python dump_model_scores.py
```

You can generate a summary of the model results as a tsv file:

```bash
python summarize_model_results.py
```

### Statistical analysis

To run the statistical analysis:

```bash
Makefile.analysis
```

### Figures

To generate all the figures at once:

```bash
make -f Makefile.figures
```

To generate specific figure types, you can run individual scripts, for example,
to generate the violn plots:

```bash
Rscript generate_violin_plots.R
```

To generate time-dependent cumulative/dynamic AUC plots:

```bash
python generate_td_auc_plots.py
```

To generate ROC and PR curve plots:

```bash
python generate_roc_pr_plots.py
```
