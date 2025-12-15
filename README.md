# NOTE: THIS REPRODUCIBLE CODE IS CURRENTLY UNFINISHED WHILE WE TEST THE PIPELINE ON OUR LOCAL CLUSTER.

# Introduction

This repository contains the code used to perform the analyses in the paper "Mind the gap: characterizing bias due to population mismatch in two-sample Mendelian randomization". It can be divided into three sections: 

1. Performing two-sample MR (2SMR) across 128 exposures and 88 outcomes.

2. Downstream processing of 2SMR results, such as Steiger filtering and estimate rescaling

3. Downstream analyses on 2SMR results such as calculating concordance z-scores, calculating shrinkage coefficients using simulation extrapolation (SIMEX), and the replication analysis assessing the tradeoff between bias and power in Biobank Japan.

# Pre-requisite Software

You will need a conda installation and R. We used R 4.5.2 for this analysis. We have tested this repository in Linux only. If you want to run the whole analysis you will also need a compute cluster. Overall, this empirical analysis will involve downloading ~140GB of GWAS summary statistics, ~10GB of LD reference panel data, and will produce an additional ~2.5 GB of output files.

To set up a

# Step 1: Performing Two-Sample MR

The packages for performing MR are located in the `TwoSampleMR`, `GRAPPLE`, and `MRBEE` packages. All of these packages are available on Github at the following links:

1. TwoSampleMR: https://github.com/MRCIEU/TwoSampleMR 
2. GRAPPLE: https://github.com/jingshuw/GRAPPLE 
3. MRBEE: https://github.com/noahlorinczcomi/MRBEE 

Other necessary packages for running the Snakemake pipeline include: VariantAnnotation (https://www.bioconductor.org/packages/release/bioc/html/VariantAnnotation.html), gwasvcf (https://github.com/MRCIEU/gwasvcf), dplyr (CRAN), data.table (CRAN), ieugwasr (https://github.com/MRCIEU/ieugwasr), and GFA (https://github.com/jean997/GFA). 

## Step 1a) Downloading summary statistics

The GWAS summary statistics used in our paper come from either OpenGWAS (which will typically be VCF files), the NHGRI-EBI GWAS catalog (which are usually non-VCF), or Release 10 of FinnGen (which are also usually non-VCF). All of these summary statistics are publicly available from their respective websites, but we provide a Zenodo record of these summary statistics for ease of use, located [here](https://zenodo.org/uploads/17584746).

This Zenodo record is divided into two files: `exposure_sumstats.tar.gz` (\~80GB) consisting of summary statistics from 128 GWAS studies across 26 exposure traits, and `outcome_sumstats.tar.gz` (\~50GB) consisting of summary statistics from 88 GWAS studies across 21 outcome traits. Once these files are downloaded and placed into your working directory, unzipping these .tar.gz files using the commands `tar -xvzf exposure_sumstats.tar.gz` and `tar -xvzf outcome_sumstats.tar.gz` will yield an "exposures/" and "outcomes/" directory, respectively, with subdirectories containing GWAS summary statistics for studies of specific traits.

## Step 1b) Downloading LD reference panels

We used reference panels from 1000 Genomes, divided into five continental superpopulations (AFR, AMR, EAS, EUR, SAS). These reference panels are available for FTP download. To obtain these summary statistics, enter the command ```wget http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz```. Then, create a folder named "ld_reference" in your working directory, and in the directory, unpack the resulting `.tgz` file using the command `tar -xvzf 1kg.v3.tgz -C ld_reference/`.

## Step 1c) Running GRAPPLE at 1e-5 with Snakemake

We use a slightly modified version of Dr. Jean Morrison's Snakemake pipeline for GFA/MVMR/NESMR located [in this repository](https://github.com/jean997/gfa_pipeline), with a focus on the multivariable MR (MVMR) portion.

The necessary components for running are provided in the following directories:

1. CSV files for each pair of exposure and outcome GWAS summary statistics, located in the `snakemake_csv_files/trait_pair_csvs/` directory.
2. The LD reference panels downloaded and unpacked from Step 1b). 
3. YAML files for each of 128 exposure studies, located in the `all_trait_yaml/` directory. Each of these YAML files instructs the snakemake pipeline to compute up to 88 MR estimates, using the specified study as an exposure across all outcome studies.

To run the snakemake pipeline for a specific exposure study:

1. Navigate to the /`snakemake-pipeline/` directory.
2. Replace the contents of `snakemake-pipeline/config_mvmr.yaml` with the corresponding YAML file in `all_trait_yaml/` (for example: `snakemake-pipeline/all_trait_yaml/bbj-a-12_mvmr.yaml`).
3. Run the command `./run-snakemake-mvmr.sh`.

The results of this snakemake output are saved into an RDS file into a data directory, `snakemake_outputs/mvmr_data/` and an output directory, `snakemake_outputs/mvmr_results`. If the snakemake pipeline is run successfully, the data directory will contain estimates of residual correlations as an RDS file, while the output directory will contain the final output of the snakemake pipeline. An sample residual correlation estimate file in the data directory is provided in example1_Restimate.RDS of the snakemake-pipeline directory, consisting of the list showing the 2x2 residual correlation matrix in `R`, and the IDs of the summary statistics used in `names`. A sample final output file in the output directory is shown in example1_run.RDS of the snakemake-pipeline directory. This object contains both the summary statistics for instruments used in the 2SMR analysis in `grapple_data`, as well as the GRAPPLE object itself, runtime information, and the IDs of the summary statistics used, in `result`.


# Step 2: Downstream Processing

The estimates derived from the output from step 1 are **not** the final estimates we used in the analyses of our paper. We need to perform three additional steps: 1) Steiger filtering to address possible reverse causation in the variants we selected as instruments, 2) performing MR using other methods at different instrument selection thresholds using our Steiger-filtered list of SNPs, and 3) rescaling the estimates to a standardized scale to facilitate comparison.

## Step 2a) Steiger Filtering



## Step 2b) Additional MR analyses



## Step 2c) Estimate Rescaling

After Steiger-filtering, the estimates must be rescaled.

# Step 3: 2SMR Analyses

Now that you have Steiger-filtered and rescaled estimates for GRAPPLE at instrument selection 
