# Introduction

This repository contains the code used to perform the analyses in the paper "Mind the gap: characterizing bias due to population mismatch in two-sample Mendelian randomization". It can be divided into three sections: 

1. Performing two-sample MR (2SMR) across 128 exposures and 88 outcomes.

2. Downstream processing of 2SMR results, such as Steiger filtering and estimate rescaling

3. Downstream analyses on 2SMR results such as calculating concordance z-scores, calculating shrinkage coefficients using simulation extrapolation (SIMEX), and the replication analysis assessing the tradeoff between bias and power in Biobank Japan.

# Step 1: Performing Two-Sample MR

## Step 1a) Running GRAPPLE at 1e-5 with Snakemake

We use a slightly modified version of Dr. Jean Morrison's Snakemake pipeline for GFA/MVMR/NESMR located [in this repository](https://github.com/jean997/gfa_pipeline), with a focus on the multivariable MR (MVMR) portion. We first ran only GRAPPLE at the lenient instrument selection p-value threshold of 1e-5, then used the resulting snakemake object to run the MR analyses at a more stringent instrument selection threshold of 5e-8. Detailed instructions for running the pipeline are located in the README of the repository, but are summarized here into general steps:

1. Download summary statistics either from OpenGWAS (which will typically be VCF files) or the EBI GWAS Catalog (which will usually be non-VCF).

2. Create csv files for the trait pairs you wish to run. These CSV files should contain key information about the summary statistics downloaded in Step 1, with the outcome in the first row and the exposure in the second. For specific examples of how to format these CSVs, two examples are provided in the snakemake-pipeline directory: example_csv1.csv that uses VCF files from OpenGWAS, and example_csv2.csv that uses both a VCF file from OpenGWAS and a non-VCF file from the GWAS Catalog. Further guidance is also provided in Section 3.2 of the original repository README. Note that the summary statistics referenced in the examples themselves are not stored in this repository due to space limitations.

3. Edit the config_mvmr.yaml file. More information about specific options are described in the comments and in the original README of the repository in Section 3.3. With the exception of the input files, output directory, and `ref_path`, all options displayed in the YAML file ones we used in our analyses. `ref_path` describes the reference panels used for LD pruning. In our study, we used reference panels from 1000 Genomes corresponding to the continent (AFR, AMR, EAS, EUR, SAS) that most closely matched the exposure. These reference panels are publicly accessible and provided by the MRC Integrative Epidemiology Unit (MRC-IEU) via FTP download [here](http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz). This link was obtained from the documentation of the ieugwasr R package on Github provided by the same organization, on the “Running local LD operations” page.

Note that it is only possible to LD prune based on one reference panel at a time. This means that we can only run exposures from the same population as the reference panel. For our analysis, we created yaml files that separated out our study pairs by shared exposure studies, resulting in 128 unique yaml files that each ran MVMR based on 88 input csv files.

4. Run the pipeline using snakemake. The README of the original repository describes how to in Section 4. We also provide a shell script (run-snakemake-mvmr.sh) that allows the pipeline to be run on a computing cluster, as well as a cluster.yaml file to configure various cluster parameters (note that by default this yaml file uses SLURM).

The results of this snakemake output are saved into an RDS file into a data directory and an output directory. If the snakemake pipeline is run successfully, the data directory will contain estimates of residual correlations as an RDS file, while the output directory will contain the final output of the snakemake pipeline. An sample residual correlation estimate file in the data directory is provided in example1_Restimate.RDS of the snakemake-pipeline directory, consisting of the list showing the 2x2 residual correlation matrix in `R`, and the IDs of the summary statistics used in `names`. A sample final output file in the output directory is shown in example1_run.RDS of the snakemake-pipeline directory. This object contains both the summary statistics for instruments used in the 2SMR analysis in `grapple_data`, as well as the GRAPPLE object itself, runtime information, and the IDs of the summary statistics used, in `result`.

## Step 1b) Running IVW and GRAPPLE at 5e-8

# Step 2: Downstream Processing

## Steiger Filtering

## Estimate Rescaling

# Step 3: 2SMR Analyses
