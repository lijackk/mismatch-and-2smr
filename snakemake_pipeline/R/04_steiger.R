source("renv/activate.R")
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(psych)

#Defining exposure and outcome traits
exposure <- snakemake@wildcards[["exposure"]]
print(exposure)
outcome <- snakemake@wildcards[["outcome"]]
print(outcome)

#Defining input and output paths
sumstat_input_path <- snakemake@input[["sumstat_final"]]
print(sumstat_input_path)
sumstat_output_path <- snakemake@output[["steiger_sumstat"]]
print(sumstat_output_path)

#Obtaining the prevalence and trait for our outcome of interest
outcome_prevalences <- fread("outcome_prevalences.csv")
outcome_prevalences <- select(outcome_prevalences, 1:4)
names(outcome_prevalences) <- c("id", "trait", "population", "prevalence")
prevalence <- outcome_prevalences$prevalence[outcome_prevalences$id == outcome] %>% as.numeric()
outcome_trait <- outcome_prevalences$trait[outcome_prevalences$id == outcome]

#Reading in the summary statistics
sumstats_prefilter <- readRDS(sumstat_input_path)

#I noted anomalous summary statistics in ieu-b-5070 the first time around - I'm still not entirely sure what the best way to handle them is, but I do know we need to adjust for them in the same way as before
if (outcome == "ieu-b-5070") {
  print("Correcting for summary statistics for outcome study ieu-b-5070")
  sumstats_prefilter$by <- log(abs(sumstats_prefilter$by))
}

#Certain outcome studies do not have allele frequency information - we need to borrow from 1000G in these instances.
  #To facilitate this process, for each of these outcome studies I have recorded the instruments across all 128 MR estimates and recorded their allele frequencies
if (outcome %in% c("ebi-a-GCST004132", "ieu-b-18", "ieu-a-832", "ieu-a-833")) {
  
  print("OpenGWAS VCF file has no AF information, AF estimated from EUR 1000G.")
  
  AF_EUR_1KG <- fread(paste0("special_afs/EUR_B38_freq_", outcome, ".frq"))
  
  AF_data <- data.frame(SNP = sumstats_prefilter$SNP)
  AF_data <- left_join(AF_data, AF_EUR_1KG, by = "SNP")
  sumstats_prefilter$af_y <- AF_data$MAF
  #AF_data <- AF_data[AF_data$MAF > 0,] #in case monomorphic alleles become a problem
}

#Collecting sample size information for our exposure and outcome
exposure_info <- fread("study_info/exposure_info.csv")
outcome_info <- fread("study_info/outcome_info.csv")

names(exposure_info) <- c("trait", "name", "pop.general", "pop.detailed", "descriptive.name", "pop.majority", "n", "units", "mean", "sd", "units.standard", "mean.standardized", "sd.standard", "transformation", "scalable", "scaling.factor")
names(outcome_info) <- c("trait", "name", "pop.general", "pop.detailed", "descriptive.name", "pop.majority", "ncase", "ncontrol", "n", "ukbb.rescale")

exposure_info$n <- as.numeric(gsub(",", "", exposure_info$n))
outcome_info$ncase <- as.numeric(gsub(",", "", outcome_info$ncase))
outcome_info$ncontrol <- as.numeric(gsub(",", "", outcome_info$ncontrol))
outcome_info$n <- as.numeric(gsub(",", "", outcome_info$n))

#Then, we run Steiger filtering, keeping only IVs where rsq.exposure is not exceeded by rsq.outcome.

sumstats_postfilter <- sumstats_prefilter
sumstats_postfilter <- sumstats_postfilter[sumstats_postfilter$byse < 100,] #stability issues.

#Sample size information
sumstats_postfilter$exposure.n <- exposure_info$n[exposure_info$name == exposure]
sumstats_postfilter$outcome.ncase <- outcome_info$ncase[outcome_info$name == outcome]
sumstats_postfilter$outcome.ncontrol <- outcome_info$ncontrol[outcome_info$name == outcome]
sumstats_postfilter$rsq.exposure <- 0
sumstats_postfilter$rsq.outcome <- 0

#Calculating rsq.exposure (a one-liner)
sumstats_postfilter$rsq.exposure <- (sumstats_postfilter$bx/sumstats_postfilter$bxse)^2 / ((sumstats_postfilter$bx/sumstats_postfilter$bxse)^2 - 2 + sumstats_postfilter$exposure.n)

#Calculating rsq.outcome
if (outcome_trait == "Total bilirubin") { #our only continuous outcome
  sumstats_postfilter$outcome.neff <- sumstats_postfilter$outcome.ncase
  sumstats_postfilter$rsq.outcome <- (sumstats_postfilter$by/sumstats_postfilter$byse)^2 / ((sumstats_postfilter$by/sumstats_postfilter$byse)^2 - 2 + sumstats_postfilter$outcome.neff)
} else {
  sumstats_postfilter <- sumstats_postfilter[!is.na(sumstats_postfilter$af_y),] #we can only calculate rsq.outcome for SNPs for which we have allele frequency information
  sumstats_postfilter$outcome.neff <- 2/(1/sumstats_postfilter$outcome.ncase + 1/sumstats_postfilter$outcome.ncontrol)
  sumstats_postfilter$rsq.outcome <- TwoSampleMR::get_r_from_lor(lor = sumstats_postfilter$by, af = sumstats_postfilter$af_y,
                                                                      ncase = sumstats_postfilter$outcome.ncase, ncontrol = sumstats_postfilter$outcome.ncontrol, prevalence = prevalence)^2
}

sumstats_postfilter$keep <- sumstats_postfilter$rsq.exposure > sumstats_postfilter$rsq.outcome

saveRDS(sumstats_postfilter, file = sumstat_output_path)