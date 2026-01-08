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

This Zenodo record is divided into two files: `exposure_sumstats.tar.gz` (\~80GB) consisting of summary statistics from 128 GWAS studies across 26 exposure traits, and `outcome_sumstats.tar.gz` (\~50GB) consisting of summary statistics from 88 GWAS studies across 21 outcome traits. Once these files are downloaded and placed into your working directory, unzipping these .tar.gz files using the commands `tar -xzvf exposure_sumstats.tar.gz` and `tar -xzvf outcome_sumstats.tar.gz` will yield an "exposures/" and "outcomes/" directory, respectively, with subdirectories containing GWAS summary statistics for studies of specific traits.

## Step 1b) Downloading LD reference panels

We used reference panels from 1000 Genomes, divided into five continental superpopulations (AFR, AMR, EAS, EUR, SAS). These reference panels are available for FTP download. To obtain these summary statistics, enter the command ```wget http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz```. Then, unpack the resulting `.tgz` file into the provided `ld_reference` directory using the command `tar -xvzf 1kg.v3.tgz -C ld_reference/`.

# Step 2: Performing two-sample MR

After setting up the data, cd into the `snakemake_pipeline` directory. The Snakefile and R scripts in this directory perform two-sample MR between the 128 exposure studies described in `snakemake_csv_files/all_exposures.csv` and the 88 outcome studies in `snakemake_csv_files/all_outcomes.csv`. Because Snakemake struggles with DAG creation when all exposures and outcomes are run together, the default behavior of the snakemake pipeline is to perform two-sample MR between all exposure studies and a single outcome study, specified in `snakemake_csv_files/current_outcome.csv`. The `run_input.R` script automatically updates `snakemake_csv_files/current_outcome.csv` across all outcomes, and submits the jobs implied by the snakemake pipeline to the cluster using the `run-snakemake.sh` command once the outcome is updated. For each individual outcome, 6,017 jobs will be sent to the job scheduler.

To run the pipeline, start an R session in the `snakemake_pipeline` directory, and paste the contents of `run_input.R` into the console. You can also use commands such as `Rscript`. Please note that by default, `run_input.R` assumes that this analysis pipeline is being performed on a computing cluster with the **Slurm job scheduler**. During our tests on our own computing cluster, and using the cluster parameters specified by both `run-snakemake.sh` and `cluster.yaml`, it took approximately **5 days** to run the analysis from start to completion, with MR analyses for individual outcomes being completed in approximately **1-2 hours each**. If your cluster uses a different type of job scheduler, you must edit the `run-snakemake.sh` script accordingly. Due to the length of this analysis, it is highly recommended to run this pipeline using either `tmux` or `screen` so that it can run in a separate background window.

# Step 3: Reproducing figures

The `reproduce_figures.R` script in the `snakemake_pipeline` directory contains all code needed to reproduce the figures and tables in both the paper and its supplement. It will output these results into the `reproducible_figures` folder in the `snakemake_pipeline` directory. We also provide expected outputs for these figures and tables in the `expected_figure_outputs` directory.
