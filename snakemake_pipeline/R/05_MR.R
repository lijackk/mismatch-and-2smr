## Running MR-IVW, GRAPPLE, and BEE on the Steiger object
  #We rescale in step 5.
source("renv/activate.R")
library(MRBEE)
library(GRAPPLE)
library(TwoSampleMR)
library(dplyr)
library(data.table)

sumstats_postfilter_path <- snakemake@input[["steiger_sumstat"]]
print(sumstats_postfilter_path)
Restimate_path <- snakemake@input[["Restimate"]]
print(Restimate_path)

mr_object_path <- snakemake@output[["mr_results"]]
print(mr_object_path)

sumstats_postfilter <- readRDS(sumstats_postfilter_path)
Restimate <- readRDS(Restimate_path)

sumstats_postfilter1e5 <- sumstats_postfilter[sumstats_postfilter$keep == TRUE,]
sumstats_postfilter5e8 <- sumstats_postfilter[sumstats_postfilter$keep == TRUE & sumstats_postfilter$selection_pvals < 5e-8,]

#We only run MR on trait pairs that have instruments left after Steiger filtering, starting with 1e-5
if (dim(sumstats_postfilter1e5)[1] > 1) {
  #=GRAPPLE=
  grapple_data1e5 <- data.frame(SNP = sumstats_postfilter1e5$SNP,
                                gamma_exp1 = sumstats_postfilter1e5$bx,
                                gamma_out = sumstats_postfilter1e5$by,
                                se_exp1 = sumstats_postfilter1e5$bxse,
                                se_out = sumstats_postfilter1e5$byse,
                                selection_pvals = sumstats_postfilter1e5$selection_pvals)
  
  #We only run diagnostics if there are more than 7 instruments
  mr_grapple_object1e5 <- grappleRobustEst(data = grapple_data1e5,
                                           p.thres = 1e-5, 
                                           plot.it = FALSE,
                                           cor.mat = Restimate$R,
                                           diagnosis = dim(grapple_data1e5)[1] > 7,
                                           loss.function = "tukey",
                                           niter = 20)
   
  #=BEE=
  mr_bee_object1e5 <- MRBEE.IMRP.UV(by = sumstats_postfilter1e5$by,
                             bx = sumstats_postfilter1e5$bx,
                             byse = sumstats_postfilter1e5$byse,
                             bxse = sumstats_postfilter1e5$bxse,
                             Rxy = Restimate$R, pv.thres = -0.02, FDR = T, var.est = "robust")
                             
  #==We only run this downstream chunk if there are at least 2 instruments with p < 5e-8==
  if (dim(sumstats_postfilter5e8)[1] > 1) {
      #=GRAPPLE
    grapple_data5e8 <- data.frame(SNP = sumstats_postfilter5e8$SNP,
                                gamma_exp1 = sumstats_postfilter5e8$bx,
                                gamma_out = sumstats_postfilter5e8$by,
                                se_exp1 = sumstats_postfilter5e8$bxse,
                                se_out = sumstats_postfilter5e8$byse,
                                selection_pvals = sumstats_postfilter5e8$selection_pvals)
                                
    #We only run diagnostics if there are more than 7 instruments
    mr_grapple_object5e8 <- grappleRobustEst(data = grapple_data5e8,
                                         p.thres = 5e-8, 
                                         plot.it = FALSE,
                                         cor.mat = Restimate$R,
                                         diagnosis = dim(grapple_data5e8)[1] > 7,
                                         loss.function = "tukey",
                                         niter = 20)
    #=BEE=
    mr_bee_object5e8 <- MRBEE.IMRP.UV(by = sumstats_postfilter5e8$by,
                               bx = sumstats_postfilter5e8$bx,
                               byse = sumstats_postfilter5e8$byse,
                               bxse = sumstats_postfilter5e8$bxse,
                               Rxy = Restimate$R, pv.thres = -0.02, FDR = T, var.est = "robust")                                           
    #=IVW=
    mr_ivw_object <- mr_ivw(b_exp = sumstats_postfilter5e8$bx,
                          b_out = sumstats_postfilter5e8$by,
                          se_exp = sumstats_postfilter5e8$bxse,
                          se_out = sumstats_postfilter5e8$byse)                                         
  } else {
    mr_grapple_object5e8 <- NA
    mr_bee_object5e8 <- NA
    mr_ivw_object <- NA
  }   
} else {
  mr_grapple_object1e5 <- NA
  mr_bee_object1e5 <- NA
  mr_grapple_object5e8 <- NA
  mr_bee_object5e8 <- NA
  mr_ivw_object <- NA  
}

#Running GRAPPLE without the filtering step for supplemental data and figures

if (dim(sumstats_postfilter)[1] > 1) {
  grapple_data1e5 <- data.frame(SNP = sumstats_postfilter$SNP,
                                    gamma_exp1 = sumstats_postfilter$bx,
                                    gamma_out = sumstats_postfilter$by,
                                    se_exp1 = sumstats_postfilter$bxse,
                                    se_out = sumstats_postfilter$byse,
                                    selection_pvals = sumstats_postfilter$selection_pvals) 
  
  grapple_data5e8 <- grapple_data1e5[grapple_data1e5$selection_pvals < 5e-8,]
  
  #We only run diagnostics if there are more than 7 instruments
  mr_grapple_object1e5.unfiltered <- grappleRobustEst(data = grapple_data1e5,
                                           p.thres = 1e-5, 
                                           plot.it = FALSE,
                                           cor.mat = Restimate$R,
                                           diagnosis = dim(grapple_data1e5)[1] > 7,
                                           loss.function = "tukey",
                                           niter = 20)
                             
  #==We only run this downstream chunk if there are at least 2 instruments with p < 5e-8==
  if (dim(grapple_data5e8)[1] > 1) {
                                
    #We only run diagnostics if there are more than 7 instruments
    mr_grapple_object5e8.unfiltered <- grappleRobustEst(data = grapple_data5e8,
                                         p.thres = 5e-8, 
                                         plot.it = FALSE,
                                         cor.mat = Restimate$R,
                                         diagnosis = dim(grapple_data5e8)[1] > 7,
                                         loss.function = "tukey",
                                         niter = 20)                                  
  } else {
    mr_grapple_object5e8.unfiltered <- NA
  }   
} else {
  mr_grapple_object1e5.unfiltered <- NA
  mr_grapple_object5e8.unfiltered <- NA
}


#Finally, we need to save the rescaling factors for the estimate, if it exists
exposure_info <- fread("study_info/exposure_info.csv")
outcome_info <- fread("study_info/outcome_info.csv")

exposure <- snakemake@wildcards[["exposure"]]
outcome <- snakemake@wildcards[["outcome"]]

names(exposure_info) <- c("trait", "name", "pop.general", "pop.detailed", "descriptive.name", "pop.majority", "n", "units", "mean", "sd", "units.standard", "mean.standardized", "sd.standard", "transformation", "scalable", "scaling.factor")
names(outcome_info) <- c("trait", "name", "pop.general", "pop.detailed", "descriptive.name", "pop.majority", "ncase", "ncontrol", "n", "ukbb.rescale")

exposure_info$n <- as.numeric(gsub(",", "", exposure_info$n))
outcome_info$ncase <- as.numeric(gsub(",", "", outcome_info$ncase))
outcome_info$ncontrol <- as.numeric(gsub(",", "", outcome_info$ncontrol))
outcome_info$n <- as.numeric(gsub(",", "", outcome_info$n))

mu <- outcome_info$ncase[outcome_info$name == outcome]/outcome_info$n[outcome_info$name == outcome]

#We divide by both the exposure and outcome rescaling factor to get our final estimate
exposure.rescale <- exposure_info$scaling.factor[exposure_info$name == exposure]
outcome.rescale <- ifelse(outcome_info$ukbb.rescale[outcome_info$name == outcome] == "Yes", mu * (1 - mu), 1)

saveRDS(list(dat1e5 = sumstats_postfilter1e5, dat5e8 = sumstats_postfilter5e8,
             mr_grapple_object1e5 = mr_grapple_object1e5, mr_bee_object1e5 = mr_bee_object1e5,
             mr_grapple_object5e8 = mr_grapple_object5e8, mr_bee_object5e8 = mr_bee_object5e8, mr_ivw_object = mr_ivw_object,
             mr_grapple_object1e5.unfiltered = mr_grapple_object1e5.unfiltered, mr_grapple_object5e8.unfiltered = mr_grapple_object5e8.unfiltered,
             exposure.rescale = exposure.rescale, outcome.rescale = outcome.rescale),
        file = mr_object_path)