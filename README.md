# NOTE: THIS REPRODUCIBLE CODE IS CURRENTLY UNFINISHED

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

This Zenodo record is divided into two files: `exposure_sumstats.tar.gz` (\~80GB) consisting of summary statistics from 128 GWAS studies across 26 exposure traits, and `outcome_sumstats.tar.gz` (\~50GB) consisting of summary statistics from 88 GWAS studies across 21 outcome traits. Once these files are downloaded and placed into your working directory, unzipping these .tar.gz files using the commands `tar -xvzf exposure_sumstats.tar.gz` and `tar -xvzf outcome_sumstats.tar.gz` will yield an "exposures/" and "outcomes/" directory, respectively, with subdirectories containing GWAS summary statistics for studies of specific traits.

## Step 1b) Downloading LD reference panels

We used reference panels from 1000 Genomes, divided into five continental superpopulations (AFR, AMR, EAS, EUR, SAS). These reference panels are available for FTP download. To obtain these summary statistics, enter the command ```wget http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz```. Then, create a folder named "ld_reference" in your working directory, and in the directory, unpack the resulting `.tgz` file using the command `tar -xvzf 1kg.v3.tgz -C ld_reference/`.

## Step 2: Performing two-sample MR

The snakemake pipeline

# Step 3: Reproducing figures
