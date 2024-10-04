# tcga-microbiome-prediction

Software: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5838055.svg)](https://doi.org/10.5281/zenodo.5838055)
Dataset: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5221525.svg)](https://doi.org/10.5281/zenodo.5221525)

[Hermida, L.C., Gertz, E.M. & Ruppin, E. Predicting cancer prognosis and drug
response from the tumor microbiome. Nat Commun 13, 2896 (2022).](https://doi.org/10.1038/s41467-022-30512-3)

To reproduce the work associated with this project, please follow the steps
below in order.

### Prerequisites

The project was developed under GNU Linux and MacOS and assumes the
use of a Unix command line shell. Both Linux and MacOS provide a
command line shell by default. One must also install developer tools
for the system, at a minimum git, make and a C/C++ compiler. Other
needed tools will be installed by the instructions below.

### Installation

Install and set up
[Miniforge3](https://github.com/conda-forge/miniforge#miniforge3)

Miniforge3 is designed to be small, and the installation will take 1-5 minutes on a
typical computer, depending on the internet connection.

Miniforge3 supplies mamba and conda, respectively, a tools for managing project
software dependencies and setting up a reproducible environment in which to run code.

`mamba` is a reimplementation and drop-in replacement for `conda`, written in
C++ and offering much faster performance and more reliable dependency solving.

To obtain the source of the project and create a conda environment
with the tools needed to run the project, execute the following below (if
using Miniforge3 replace `mamba` with `conda`):

```bash
git clone --recurse-submodules https://github.com/ruppinlab/tcga-microbiome-prediction.git
cd tcga-microbiome-prediction
mamba env create -f envs/tcga-microbiome-prediction.yml
mamba activate tcga-microbiome-prediction
```

The time to complete this step is dependent on your internet
connection, but with a typical computer and internet connection takes 1-2
minutes. The newly created conda environment installs several
software packages, listed with their version numbers in
`envs/tcga-microbiome-prediction.yml`. In particular, it contains:

- python 3.10.14
- R 4.3.3

These packages are only visible within the active conda environment and
`mamba activate tcga-microbiome-prediction` only applies to the current command
line shell where it was activated.

### Data preprocessing

Process the Kraken2 + Bracken microbial abundance data
(originally generated from our pipeline
[tcga-wgs-kraken-microbial-quant](https://github.com/hermidalc/tcga-wgs-kraken-microbial-quant)
) as well as relevant case and sample metadata data from the NCI Genomic Data Commons
(GDC):

```bash
Rscript process_k2b_data_gdc_meta.R
```

Get gene annotation data from GENCODE and Ensembl for GENCODE v36 and Ensembl v102:

```bash
Rscript get_gtf_ensg_annots.R
```

Get latest TCGA patient drug response data from the NCI GDC and map to standard
drug names:

```bash
Rscript get_drug_response_data.R
```

Process TCGA survival and drug response phenotypic data:

```bash
Rscript process_surv_resp_pdata.R
```

Download NCI GDC gene expression data and create microbial abundance, gene
expression, and combined datasets in a format appropriate for the ML code
(ExpressionSet objects). This step takes ~3GB space:

```bash
Rscript create_k2b_esets.R
```

The most significant time in the data preprocessing step is typically
when running the `create_k2b_esets.R` script. The actual time required will
depend not only on the speed of the internet connection, but on how
busy the NCI GDC servers are. With typical internet connection speeds it
takes ~1 hour.

### Clinical covariate ML models

Create the clinical covariate-only ML models -- models containing age at
diagnosis, gender and tumor stage. Create these and save the results:

```bash
python run_surv_clinical_models.py
python run_resp_clinical_models.py
```

### Microbial abundance, gene expression, and combination data type ML models

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

The modeling results will be in `results/models`.

### Model results

Extract and save model scores in formats natural for python and GNU R
(pandas and R dataframes):

```bash
python dump_model_scores.py
```

You can generate a summary of the model results as a TSV file:

```bash
python summarize_model_results.py
```

### Statistical analysis

To run the statistical analysis:

```bash
make -f Makefile.analysis
```

The output of the analysis will be in `results/analysis`. See the
README therein for a description of the output files.

### Figures

To generate all the figures at once:

```bash
make -f Makefile.figures
```

The figures will be created as subdirectories under the directory `figures`.
Please note that some figure type scripts take ~30 minutes each to complete.

To generate specific figure types, you can run individual scripts, for example,
to generate the violin plots:

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

To generate permutation test histogram plots:

```bash
python generate_perm_hist_plots.py
```
