# NOTE: THIS REPRODUCIBLE CODE IS A WORK IN PROGRESS

# Introduction

This repository contains the code used to perform the analyses in the paper "Mind the gap: characterizing bias due to population mismatch in two-sample Mendelian randomization". It can be divided into three sections: 

1. Downloading required data (MR summary statistics and LD reference panels)

2. Performing two-sample MR (2SMR) across 128 exposures and 88 outcomes.

3. Downstream analyses on 2SMR results, replicating results and figures from the paper.

# Step 1: Pre-requisite Software

You will need a conda installation and R. We used R 4.5.2 for this analysis. We have tested this repository in Linux only. If you want to run the whole analysis you will also need a compute cluster. Overall, this empirical analysis will involve downloading ~140GB of GWAS summary statistics, ~10GB of LD reference panel data, and will produce an additional ~4.2 GB of output files.

All necessary packages and the versions used in our analyses are standardized using the `renv` folder in the `snakemake-pipeline` directory.

## Step 1a) Downloading summary statistics

The GWAS summary statistics used in our paper come from either OpenGWAS (which will typically be VCF files), the NHGRI-EBI GWAS catalog (which are usually non-VCF), or Release 10 of FinnGen (which are also usually non-VCF). All of these summary statistics are publicly available from their respective websites, but we provide a Zenodo record of these summary statistics for ease of use, located [here](https://zenodo.org/uploads/17584746).

This Zenodo record is divided into two files: `exposure_sumstats.tar.gz` (\~80GB) consisting of summary statistics from 128 GWAS studies across 26 exposure traits, and `outcome_sumstats.tar.gz` (\~50GB) consisting of summary statistics from 88 GWAS studies across 21 outcome traits. Once these files are downloaded and placed into your working directory, unzipping these .tar.gz files using the commands

```
tar -xzvf exposure_sumstats.tar.gz
tar -xzvf outcome_sumstats.tar.gz
```

will yield an `exposures/` and `outcomes/` directory, respectively, with subdirectories containing GWAS summary statistics for studies of specific traits.

## Step 1b) Downloading LD reference panels

We used reference panels from 1000 Genomes, divided into five continental superpopulations (AFR, AMR, EAS, EUR, SAS). These reference panels are available for FTP download. To obtain these summary statistics, enter the command

```
wget http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz
```

Then, unpack the resulting `.tgz` file into the provided `ld_reference` directory using the command

```
tar -xvzf 1kg.v3.tgz -C ld_reference/
```

# Step 2: Performing two-sample MR

## Step 2a: Setting up conda and R environments

After setting up the data, navigate into the `snakemake_pipeline/` directory using

```
cd snakemake_pipeline/
```

The Snakefile and R scripts in this directory perform two-sample MR between the 128 exposure studies described in `snakemake_csv_files/all_exposures.csv` and the 88 outcome studies in `snakemake_csv_files/all_outcomes.csv`. This snakemake pipeline runs in a conda environment, which needs to be built and activated. To build the environment, type

```
conda env create -f manuscript_env.yml
```

Once finished, confirm that the environment was built successfully using

```
conda env list
```

You should see `snakemake_env` in this list of environments. Activate this environment using

```
conda activate snakemake_env
```

This snakemake pipeline uses the `renv` R package to ensure consistent R package environments. Once the renv package has been installed, start an R session and type the command

```
renv::restore()
```

This will install of the necessary packages for our analysis pipeline.

## Step 2b: Running the snakemake pipeline

During our tests on our own computing cluster, and using the cluster parameters specified by both `run-snakemake.sh` and `cluster.yaml`, it took approximately **5 days** to run the analysis from start to completion. Due to the length of this analysis, it is highly recommended to run this pipeline using either `tmux` or `screen` so that it can run in a separate background window.

Once in this background window, the entire pipeline can be run using

```
./run-snakemake.sh
```

Please note that by default, `run-snakemake.sh` assumes that this analysis pipeline is being performed on a computing cluster with the **Slurm job scheduler**. If your cluster uses a different type of job scheduler, you must edit the `run-snakemake.sh` script accordingly. Overall, this pipeline will submit __529,496__ jobs to the job scheduler.

The outputs for these analyses will be saved into `outputs`, while logs for specific jobs will be saved into a `snakemake_pipeline/logs` directory. If no troubleshooting is necessary after the pipeline is finished running, it is highly recommended to remove the `logs` folder using `rm -rf logs/`.

# Step 3: Reproducing figures

The `reproduce_figures.R` script in the `snakemake_pipeline` directory contains all code needed to reproduce the figures and tables in both the paper and its supplement. It will output these results into the `reproducible_figures` folder in the `snakemake_pipeline` directory. We also provide expected outputs for these figures and tables in the `expected_figure_outputs` directory. This script can be run by typing the following command into the command line:

```
Rscript reproduce_figures.R
```
