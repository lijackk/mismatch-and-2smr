#A script that reproduces the figures and key results from our analyses
source("renv/activate.R")
library(data.table)
library(dplyr)
library(purrr)
library(ggplot2)

##===Preparatory analyses===

#==Reading in the MR results and creating a table==

#Reading estimates
new_mr_paths <- list.files("../outputs/MRest", full.names = TRUE)
new_mr <- list()
for (i in 1:length(new_mr_paths)) {
  new_mr[[i]] <- readRDS(new_mr_paths[i])
  if (i %% 100 == 0) print(i)
}
new_mr_studypairs <- gsub("(../outputs/MRest/)(.*)(\\.MRestimate_list.RDS)", "\\2", new_mr_paths)
names(new_mr) <- new_mr_studypairs

#Assembling table
new.grapple1e5.beta <- vector()
new.bee1e5.beta <- vector()

new.grapple5e8.beta <- vector()
new.bee5e8.beta <- vector()
new.ivw.beta <- vector()

new.grapple1e5.se <- vector()
new.bee1e5.se <- vector()

new.grapple5e8.se <- vector()
new.bee5e8.se <- vector()
new.ivw.se <- vector()

for (i in 1:length(new_mr_studypairs)) {
  current.studypair <- new_mr_studypairs[i]

  #==Extracting MR objects==
  new.grapple1e5.object <- new_mr[[which(names(new_mr) == current.studypair)]]$mr_grapple_object1e5
  new.bee1e5.object <- new_mr[[which(names(new_mr) == current.studypair)]]$mr_bee_object1e5

  new.grapple5e8.object <- new_mr[[which(names(new_mr) == current.studypair)]]$mr_grapple_object5e8
  new.bee5e8.object <- new_mr[[which(names(new_mr) == current.studypair)]]$mr_bee_object5e8
  new.ivw.object <- new_mr[[which(names(new_mr) == current.studypair)]]$mr_ivw_object

  #==Extracting beta hat==
  new.grapple1e5.beta[i] <- ifelse(identical(new.grapple1e5.object, NA), NA, new.grapple1e5.object$beta.hat)
  new.bee1e5.beta[i] <- ifelse(identical(new.bee1e5.object, NA), NA, new.bee1e5.object$theta)

  new.grapple5e8.beta[i] <- ifelse(identical(new.grapple5e8.object, NA), NA, new.grapple5e8.object$beta.hat)
  new.bee5e8.beta[i] <- ifelse(identical(new.bee5e8.object, NA), NA, new.bee5e8.object$theta)
  new.ivw.beta[i] <- ifelse(identical(new.ivw.object, NA), NA, new.ivw.object$b)

  #==Extracting standard errors==
  new.grapple1e5.se[i] <- ifelse(identical(new.grapple1e5.object, NA), NA, new.grapple1e5.object$beta.var %>% sqrt())
  new.bee1e5.se[i] <- ifelse(identical(new.bee1e5.object, NA), NA, new.bee1e5.object$vartheta %>% sqrt())

  new.grapple5e8.se[i] <- ifelse(identical(new.grapple5e8.object, NA), NA, new.grapple5e8.object$beta.var %>% sqrt())
  new.bee5e8.se[i] <- ifelse(identical(new.bee5e8.object, NA), NA, new.bee5e8.object$vartheta %>% sqrt())
  new.ivw.se[i] <- ifelse(identical(new.ivw.object, NA), NA, new.ivw.object$se)

  if (i %% 1000 == 0) print(i)
}

mr_table <- data.frame(study_pair = new_mr_studypairs,
                       exposure_id = gsub("(.*)__(.*)", "\\1", new_mr_studypairs),
                       outcome_id = gsub("(.*)__(.*)", "\\2", new_mr_studypairs),
                       grapple1e5.beta = new.grapple1e5.beta, grapple1e5.se = new.grapple1e5.se,
                       bee1e5.beta = new.bee1e5.beta, bee1e5.se = new.bee1e5.se,
                       grapple5e8.beta = new.grapple5e8.beta, grapple5e8.se = new.grapple5e8.se,
                       bee5e8.beta = new.bee5e8.beta, bee5e8.se = new.bee5e8.se,
                       ivw.beta = new.ivw.beta, ivw.se = new.ivw.se)

mr_table$exposure.rescale <- lapply(new_mr, function(X){X$exposure.rescale}) %>% unlist()
mr_table$outcome.rescale <- lapply(new_mr, function(X){X$outcome.rescale}) %>% unlist()

#==Obtaining rescaled estimates==

mr_table$grapple1e5.beta.rescale <- mr_table$grapple1e5.beta/mr_table$exposure.rescale/mr_table$outcome.rescale
mr_table$grapple1e5.se.rescale <- mr_table$grapple1e5.se/mr_table$exposure.rescale/mr_table$outcome.rescale

mr_table$bee1e5.beta.rescale <- mr_table$bee1e5.beta/mr_table$exposure.rescale/mr_table$outcome.rescale
mr_table$bee1e5.se.rescale <- mr_table$bee1e5.se/mr_table$exposure.rescale/mr_table$outcome.rescale

mr_table$grapple5e8.beta.rescale <- mr_table$grapple5e8.beta/mr_table$exposure.rescale/mr_table$outcome.rescale
mr_table$grapple5e8.se.rescale <- mr_table$grapple5e8.se/mr_table$exposure.rescale/mr_table$outcome.rescale

mr_table$bee5e8.beta.rescale <- mr_table$bee5e8.beta/mr_table$exposure.rescale/mr_table$outcome.rescale
mr_table$bee5e8.se.rescale <- mr_table$bee5e8.se/mr_table$exposure.rescale/mr_table$outcome.rescale

mr_table$ivw.beta.rescale <- mr_table$ivw.beta/mr_table$exposure.rescale/mr_table$outcome.rescale
mr_table$ivw.se.rescale <- mr_table$ivw.se/mr_table$exposure.rescale/mr_table$outcome.rescale

#==Applying ETOS groups to the MR table dataset==
to_join <- fread("study_info/estimate_groups.csv")

mr_table.ETOS <- inner_join(mr_table, to_join, by = c("exposure_id", "outcome_id"))
mr_table.ETOS <- arrange(mr_table.ETOS, ETOS)
ETOS.with.exact <- unique(mr_table.ETOS$ETOS[mr_table.ETOS$is.bestmatch & mr_table.ETOS$match.type == "Exact Match"])
mr_table.exactETOS <- mr_table.ETOS[mr_table.ETOS$ETOS %in% ETOS.with.exact,]

#==Calculating z-scores of concordance==
mr_table.exactETOS$zconcord.grapple1e5 <- NA
mr_table.exactETOS$zconcord.bee1e5 <- NA
mr_table.exactETOS$zconcord.grapple5e8 <- NA
mr_table.exactETOS$zconcord.bee5e8 <- NA
mr_table.exactETOS$zconcord.ivw <- NA

for (i in 1:dim(mr_table.exactETOS)[1]) {
  if (mr_table.exactETOS$is.bestmatch[i]) { #No need to compare the reference estimate to itself!
    next
  }

  #Identifying the reference estimate pair
  ETOS <- mr_table.exactETOS$ETOS[i]
  ETOS.outcome <- mr_table.exactETOS$outcome_id[i]
  ETOS.reference <- mr_table.exactETOS$exposure_id[mr_table.exactETOS$is.bestmatch & mr_table.exactETOS$ETOS == ETOS & mr_table.exactETOS$outcome_id == ETOS.outcome]

  #Reference estimates and SEs for all five MR estimate types, if they exist
  grapple1e5.reference <- mr_table.exactETOS$grapple1e5.beta.rescale[mr_table.exactETOS$outcome_id == ETOS.outcome & mr_table.exactETOS$exposure_id == ETOS.reference]
  grapple1e5.reference.se <- mr_table.exactETOS$grapple1e5.se.rescale[mr_table.exactETOS$outcome_id == ETOS.outcome & mr_table.exactETOS$exposure_id == ETOS.reference]
  bee1e5.reference <- mr_table.exactETOS$bee1e5.beta.rescale[mr_table.exactETOS$outcome_id == ETOS.outcome & mr_table.exactETOS$exposure_id == ETOS.reference]
  bee1e5.reference.se <- mr_table.exactETOS$bee1e5.se.rescale[mr_table.exactETOS$outcome_id == ETOS.outcome & mr_table.exactETOS$exposure_id == ETOS.reference]
  grapple5e8.reference <- mr_table.exactETOS$grapple5e8.beta.rescale[mr_table.exactETOS$outcome_id == ETOS.outcome & mr_table.exactETOS$exposure_id == ETOS.reference]
  grapple5e8.reference.se <- mr_table.exactETOS$grapple5e8.se.rescale[mr_table.exactETOS$outcome_id == ETOS.outcome & mr_table.exactETOS$exposure_id == ETOS.reference]
  bee5e8.reference <- mr_table.exactETOS$bee5e8.beta.rescale[mr_table.exactETOS$outcome_id == ETOS.outcome & mr_table.exactETOS$exposure_id == ETOS.reference]
  bee5e8.reference.se <- mr_table.exactETOS$bee5e8.se.rescale[mr_table.exactETOS$outcome_id == ETOS.outcome & mr_table.exactETOS$exposure_id == ETOS.reference]
  ivw.reference <- mr_table.exactETOS$ivw.beta.rescale[mr_table.exactETOS$outcome_id == ETOS.outcome & mr_table.exactETOS$exposure_id == ETOS.reference]
  ivw.reference.se <- mr_table.exactETOS$ivw.se.rescale[mr_table.exactETOS$outcome_id == ETOS.outcome & mr_table.exactETOS$exposure_id == ETOS.reference]

  #The mismatching estimates and SEs for both unfiltered/Steiger
  grapple1e5.mismatch <- mr_table.exactETOS$grapple1e5.beta.rescale[i]
  grapple1e5.mismatch.se <- mr_table.exactETOS$grapple1e5.se.rescale[i]
  bee1e5.mismatch <- mr_table.exactETOS$bee1e5.beta.rescale[i]
  bee1e5.mismatch.se <- mr_table.exactETOS$bee1e5.se.rescale[i]
  grapple5e8.mismatch <- mr_table.exactETOS$grapple5e8.beta.rescale[i]
  grapple5e8.mismatch.se <- mr_table.exactETOS$grapple5e8.se.rescale[i]
  bee5e8.mismatch <- mr_table.exactETOS$bee5e8.beta.rescale[i]
  bee5e8.mismatch.se <- mr_table.exactETOS$bee5e8.se.rescale[i]
  ivw.mismatch <- mr_table.exactETOS$ivw.beta.rescale[i]
  ivw.mismatch.se <- mr_table.exactETOS$ivw.se.rescale[i]

  #Calculating concordance Z-scores
  mr_table.exactETOS$zconcord.grapple1e5[i] <- (grapple1e5.mismatch - grapple1e5.reference)/sqrt(grapple1e5.reference.se^2 + grapple1e5.mismatch.se^2)
  mr_table.exactETOS$zconcord.bee1e5[i] <- (bee1e5.mismatch - bee1e5.reference)/sqrt(bee1e5.reference.se^2 + bee1e5.mismatch.se^2)
  mr_table.exactETOS$zconcord.grapple5e8[i] <- (grapple5e8.mismatch - grapple5e8.reference)/sqrt(grapple5e8.reference.se^2 + grapple5e8.mismatch.se^2)
  mr_table.exactETOS$zconcord.bee5e8[i] <- (bee5e8.mismatch - bee5e8.reference)/sqrt(bee5e8.reference.se^2 + bee5e8.mismatch.se^2)
  mr_table.exactETOS$zconcord.ivw[i] <- (ivw.mismatch - ivw.reference)/sqrt(ivw.reference.se^2 + ivw.mismatch.se^2)
}

#==Calculating concordance z-score p-values==
mr_table.exactETOS$pconcord.grapple1e5 <- 2 * pnorm(abs(mr_table.exactETOS$zconcord.grapple1e5), lower.tail = FALSE)
mr_table.exactETOS$pconcord.grapple1e5.bh <- p.adjust(mr_table.exactETOS$pconcord.grapple1e5, method = "BH")

mr_table.exactETOS$pconcord.bee1e5 <- 2 * pnorm(abs(mr_table.exactETOS$zconcord.bee1e5), lower.tail = FALSE)
mr_table.exactETOS$pconcord.bee1e5.bh <- p.adjust(mr_table.exactETOS$pconcord.bee1e5, method = "BH")

mr_table.exactETOS$pconcord.grapple5e8 <- 2 * pnorm(abs(mr_table.exactETOS$zconcord.grapple5e8), lower.tail = FALSE)
mr_table.exactETOS$pconcord.grapple5e8.bh <- p.adjust(mr_table.exactETOS$pconcord.grapple5e8, method = "BH")

mr_table.exactETOS$pconcord.bee5e8 <- 2 * pnorm(abs(mr_table.exactETOS$zconcord.bee5e8), lower.tail = FALSE)
mr_table.exactETOS$pconcord.bee5e8.bh <- p.adjust(mr_table.exactETOS$pconcord.bee5e8, method = "BH")

mr_table.exactETOS$pconcord.ivw <- 2 * pnorm(abs(mr_table.exactETOS$zconcord.ivw), lower.tail = FALSE)
mr_table.exactETOS$pconcord.ivw.bh <- p.adjust(mr_table.exactETOS$pconcord.ivw, method = "BH")

#==Forest plots for LDL-CHD==
#=Figure 3=
ldl.chd.1e5 <- mr_table.ETOS[mr_table.ETOS$ETOS %like% "LDL cholesterol" &
                             mr_table.ETOS$outcome_id %in% c("ebi-a-GCST005194", "finngen_R10_I9_CHD", "bbj-a-159", "ieu-a-7"),] %>%
               select("exposure_id", "outcome_id", "grapple1e5.beta.rescale", "grapple1e5.se.rescale", "match.type")
ldl.chd.1e5$threshold <- "1e-5"
ldl.chd.1e5 <- left_join(ldl.chd.1e5, mr_table.exactETOS[,c(2,3,29,35)], by = c("exposure_id", "outcome_id"))
names(ldl.chd.1e5)[c(3,4,7,8)] <- c("beta.rescale", "se.rescale", "zconcord", "pconcord.bh")

ldl.chd.5e8 <- mr_table.ETOS[mr_table.ETOS$ETOS %like% "LDL cholesterol" &
                             mr_table.ETOS$outcome_id %in% c("ebi-a-GCST005194", "finngen_R10_I9_CHD", "bbj-a-159", "ieu-a-7"),] %>%
               select("exposure_id", "outcome_id", "grapple5e8.beta.rescale", "grapple5e8.se.rescale", "match.type")
ldl.chd.5e8$threshold <- "5e-8"
ldl.chd.5e8 <- left_join(ldl.chd.5e8, mr_table.exactETOS[,c(2,3,30,37)], by = c("exposure_id", "outcome_id"))
names(ldl.chd.5e8)[c(3,4,7,8)] <- c("beta.rescale", "se.rescale", "zconcord", "pconcord.bh")

ldl.chd.full <- rbind(ldl.chd.1e5, ldl.chd.5e8)

ldl.chd.full$exposure_id <- factor(ldl.chd.full$exposure_id, levels = c("bbj-a-31", "ebi-a-GCST008037", "ieu-b-110"))
levels(ldl.chd.full$exposure_id) <- c("Biobank Japan", "Majority Hispanic\nMega-analysis", "UK Biobank")
ldl.chd.full$outcome_id <- factor(ldl.chd.full$outcome_id, levels = unique(ldl.chd.full$outcome_id))
levels(ldl.chd.full$outcome_id) <- c("Biobank Japan", "UK Biobank", "FinnGen", "CARDIoGRAMplusC4D")

ldl.chd.full$match.type2 <- ifelse(ldl.chd.full$match.type == "Exact Match", "Reference Estimate", ifelse(ldl.chd.full$match.type == "Continental Mismatch", "Mismatch: Different Continent", "Mismatch: Same Continent"))

color_palette <- c("Reference Estimate" = "#0453AB", "Mismatch: Same Continent" = "#76A9E0", "Mismatch: Different Continent" = "#B10000")
shape_palette <- c("1e-5" = 16, "5e-8" = 15)

ci_pval <- 0.05
ci_width <- abs(qnorm(ci_pval/2))

ggplot(data = ldl.chd.full, aes(x = beta.rescale, y = exposure_id, col = match.type2, group = threshold, shape = threshold)) + geom_point(position = position_dodge(width=0.5), size = 3.5) + facet_wrap(~factor(paste("Target population:", outcome_id), levels = unique(paste("Target population:", outcome_id))), scales = "fixed") + geom_errorbar(aes(xmin = beta.rescale - ci_width*se.rescale, xmax = beta.rescale + ci_width*se.rescale), width = 0.25, position = position_dodge(width=0.5), data = ldl.chd.full, lwd = 1) + geom_vline(xintercept = 0, lty = 2, col = "black") + labs(y = "", col = "Population Match", x = "2SMR Estimates", shape = "IV Selection Threshold") + theme_bw(base_size = 17.5) + scale_color_manual(values = color_palette) + scale_shape_manual(values = shape_palette)
ggsave("reproducible_figures/figure3_ldlchd.png")

#==Analysis of overall shrinkage==
source("simex_source_func.R")
#=Table 1 - Concordance z-score analysis, GRAPPLE=
#Calculating "mismatch-to-reference" ratio for GRAPPLE estimates

mtr.ratio.1e5 <- vector()
reference.estimate.1e5 <- vector()
reference.se.1e5 <- vector()

mtr.ratio.5e8 <- vector()
reference.estimate.5e8 <- vector()
reference.se.5e8 <- vector()

for (i in 1:dim(mr_table.exactETOS)[1]) {
  if (mr_table.exactETOS$is.bestmatch[i]) { #mismatch/reference ratio only applies to mismatches.
    mtr.ratio.1e5[i] <- NA
    reference.estimate.1e5[i] <- NA
    reference.se.1e5[i] <- NA
    mtr.ratio.5e8[i] <- NA
    reference.estimate.5e8[i] <- NA
    reference.se.5e8[i] <- NA
    next
  }

  mismatch.estimate.1e5 <- mr_table.exactETOS$grapple1e5.beta.rescale[i]
  mismatch.estimate.5e8 <- mr_table.exactETOS$grapple5e8.beta.rescale[i]
  ETOS <- mr_table.exactETOS$ETOS[i]

  reference.estimate.1e5[i] <- mr_table.exactETOS$grapple1e5.beta.rescale[mr_table.exactETOS$is.bestmatch & mr_table.exactETOS$ETOS == ETOS]
  reference.se.1e5[i] <- mr_table.exactETOS$grapple1e5.se.rescale[mr_table.exactETOS$is.bestmatch & mr_table.exactETOS$ETOS == ETOS]
  mtr.ratio.1e5[i] <- mismatch.estimate.1e5/reference.estimate.1e5[i]

  if (!is.na(mismatch.estimate.5e8)) {
    reference.estimate.5e8[i] <- mr_table.exactETOS$grapple5e8.beta.rescale[mr_table.exactETOS$is.bestmatch & mr_table.exactETOS$ETOS == ETOS]
    reference.se.5e8[i] <- mr_table.exactETOS$grapple5e8.se.rescale[mr_table.exactETOS$is.bestmatch & mr_table.exactETOS$ETOS == ETOS]
    mtr.ratio.5e8[i] <- mismatch.estimate.5e8/reference.estimate.5e8[i]
  } else {
    reference.estimate.5e8[i] <- NA
    reference.se.5e8[i] <- NA
    mtr.ratio.5e8[i] <- NA
  }
}

#Calculating the overall shrinkage coefficients for GRAPPLE using SIMEX
set.seed(2025)
indices.1e5 <- which(!is.na(mr_table.exactETOS$zconcord.grapple1e5))
simex.grapple1e5 <- simex_regression(outcome = mr_table.exactETOS$grapple1e5.beta.rescale[indices.1e5],
                                     outcome.se = mr_table.exactETOS$grapple1e5.se.rescale[indices.1e5],
                                     predictor = reference.estimate.1e5[indices.1e5],
                                     predictor.se = reference.se.1e5[indices.1e5])

indices.5e8 <- which(!is.na(mr_table.exactETOS$zconcord.grapple5e8))
simex.grapple5e8 <- simex_regression(outcome = mr_table.exactETOS$grapple5e8.beta.rescale[indices.5e8],
                                     outcome.se = mr_table.exactETOS$grapple5e8.se.rescale[indices.5e8],
                                     predictor = reference.estimate.5e8[indices.5e8],
                                     predictor.se = reference.se.5e8[indices.5e8])


table_categories <- c("Instrument Selection Threshold",
                      "Total Number of Mismatching Estimates with Reference Estimates",
                      "Mismatching Estimates with Non-Significant Concordance z-scores",
                      "Mismatching Estimates with Significant Concordance z-scores, Same Sign as Reference, Smaller Magnitude",
                      "Mismatching Estimates with Significant Concordance z-scores, Same Sign as Reference, Larger Magnitude",
                      "Mismatching Estimates with Significant Concordance z-scores, Opposite Sign as Reference",
                      "Mismatching Estimates with Significant Concordance z-scores, Opposite Sign as Reference, Reference and Mismatch MR p-values < 0.05",
                      "SIMEX-Estimated Overall Shrinkage Coefficient, 1000 Bootstrap Resamples")

crit.z <- abs(qnorm(0.05/2))
numbers_1e5 <- c(1e-5,
                 sum(!is.na(mr_table.exactETOS$zconcord.grapple1e5)),
                 sum(mr_table.exactETOS$pconcord.grapple1e5.bh >= 0.05, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.grapple1e5.bh < 0.05 & mtr.ratio.1e5 > 0 & mtr.ratio.1e5 < 1, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.grapple1e5.bh < 0.05 & mtr.ratio.1e5 > 1, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.grapple1e5.bh < 0.05 & mtr.ratio.1e5 < 0, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.grapple1e5.bh < 0.05 & mtr.ratio.1e5 < 0 &
                     abs(mr_table.exactETOS$grapple1e5.beta.rescale/mr_table.exactETOS$grapple1e5.se.rescale) > crit.z &
                     abs(reference.estimate.1e5/reference.se.1e5) > crit.z, na.rm = TRUE),
                 paste0(round(mean(simex.grapple1e5), 3), " (", round(sd(simex.grapple1e5), 3), ")"))

numbers_5e8 <- c(5e-8,
                 sum(!is.na(mr_table.exactETOS$zconcord.grapple5e8)),
                 sum(mr_table.exactETOS$pconcord.grapple5e8.bh >= 0.05, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.grapple5e8.bh < 0.05 & mtr.ratio.5e8 > 0 & mtr.ratio.5e8 < 1, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.grapple5e8.bh < 0.05 & mtr.ratio.5e8 > 1, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.grapple5e8.bh < 0.05 & mtr.ratio.5e8 < 0, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.grapple5e8.bh < 0.05 & mtr.ratio.5e8 < 0 &
                     abs(mr_table.exactETOS$grapple5e8.beta.rescale/mr_table.exactETOS$grapple5e8.se.rescale) > crit.z &
                     abs(reference.estimate.5e8/reference.se.5e8) > crit.z, na.rm = TRUE),
                 paste0(round(mean(simex.grapple5e8), 3), " (", round(sd(simex.grapple5e8), 3), ")"))

fwrite(data.frame(table_categories, numbers_1e5, numbers_5e8), col.names = FALSE, "reproducible_figures/table1_zconcord_grapple.csv")


#=Supplemental Figure S3: Showing that a linear extrapolation in SIMEX likely works fine=
set.seed(2025)
simex.grapple1e5.nb <- simex_regression.nb(outcome = mr_table.exactETOS$grapple1e5.beta.rescale[indices.1e5],
                                     outcome.se = mr_table.exactETOS$grapple1e5.se.rescale[indices.1e5],
                                     predictor = reference.estimate.1e5[indices.1e5],
                                     predictor.se = reference.se.1e5[indices.1e5])

simex.grapple5e8.nb <- simex_regression.nb(outcome = mr_table.exactETOS$grapple5e8.beta.rescale[indices.5e8],
                                     outcome.se = mr_table.exactETOS$grapple5e8.se.rescale[indices.5e8],
                                     predictor = reference.estimate.5e8[indices.5e8],
                                     predictor.se = reference.se.5e8[indices.5e8])

png(file = "reproducible_figures/figures3_simexlinear.png", width = 1920, height = 1080)
par(mar = c(5.1, 5.1, 4.1, 2.1), mfcol = c(1,2))

plot(c(-1, seq(0, 5, by = 0.25)), c(simex.grapple1e5.nb$simex_estimate, simex.grapple1e5.nb$slope_lambda),
     main = "SIMEX Estimation of Overall Shrinkage Coefficient, GRAPPLE (1e-5)", xlab = "SIMEX Lambda", ylab = "Estimated Shrinkage Coefficient",
     cex.main = 2.2, cex.lab = 2, cex.axis = 1.8,
     pch = c(16, rep(1, 21)))
abline(a = simex.grapple1e5.nb$simex_model[1] , b = simex.grapple1e5.nb$simex_model[2], lty = 2)
abline(v = 0)

plot(c(-1, seq(0, 5, by = 0.25)), c(simex.grapple5e8.nb$simex_estimate, simex.grapple5e8.nb$slope_lambda),
     main = "SIMEX Estimation of Overall Shrinkage Coefficient, GRAPPLE (5e-8)", xlab = "SIMEX Lambda", ylab = "Estimated Shrinkage Coefficient",
     cex.main = 2.2, cex.lab = 2, cex.axis = 1.8,
     pch = c(16, rep(1, 21)))
abline(a = simex.grapple5e8.nb$simex_model[1] , b = simex.grapple5e8.nb$simex_model[2], lty = 2)
abline(v = 0)

dev.off()

#=Supplemental Table S4 - Mismatching estimates significantly larger than reference with GRAPPLE=
larger.1e5 <- which(mr_table.exactETOS$pconcord.grapple1e5.bh < 0.05 & mtr.ratio.1e5 > 1)
larger.5e8 <- which(mr_table.exactETOS$pconcord.grapple5e8.bh < 0.05 & mtr.ratio.5e8 > 1)

grapple1e5_pops <- fread("study_info/pop_id_reference.csv")
exposure_info <- fread("study_info/exposure_info.csv")
outcome_info <- fread("study_info/outcome_info.csv")

larger <- sort(c(larger.1e5, larger.5e8))
exposure_id <- mr_table.exactETOS$exposure_id[larger]
outcome_id <- mr_table.exactETOS$outcome_id[larger]

exposure_pop <- sapply(exposure_id, function(X){grapple1e5_pops$exposure_pop[grapple1e5_pops$exposure_id == X]} %>% unique()) %>% unlist() %>% unname()
outcome_pop <- sapply(outcome_id, function(X){grapple1e5_pops$outcome_pop[grapple1e5_pops$outcome_id == X]} %>% unique()) %>% unlist() %>% unname()
exposure_trait <- sapply(exposure_id, function(X){exposure_info$Trait[exposure_info$`Study ID` == X]}) %>% unlist() %>% unname()
outcome_trait <- sapply(outcome_id, function(X){outcome_info$Trait[outcome_info$`Study ID` == X]}) %>% unlist() %>% unname()

threshold <- c(5e-8, 1e-5, rep(5e-8, 3))

mismatch_estimate <- ifelse(threshold == 5e-8, mr_table.exactETOS$grapple5e8.beta.rescale[larger], mr_table.exactETOS$grapple1e5.beta.rescale[larger])
mismatch_se <- ifelse(threshold == 5e-8, mr_table.exactETOS$grapple5e8.se.rescale[larger], mr_table.exactETOS$grapple1e5.se.rescale[larger])
mismatch_p <- 2*pnorm(abs(mismatch_estimate/mismatch_se), lower.tail = FALSE)

reference_estimate <- ifelse(threshold == 5e-8, reference.estimate.5e8[larger], reference.se.1e5[larger])
reference_se <- ifelse(threshold == 5e-8, reference.se.5e8[larger], reference.se.5e8[larger])
reference_p <- 2*pnorm(abs(reference_estimate/reference_se), lower.tail = FALSE)

zconcord <- ifelse(threshold == 5e-8, mr_table.exactETOS$zconcord.grapple5e8[larger], mr_table.exactETOS$zconcord.grapple1e5[larger])

fwrite(data.frame(paste(exposure_trait, "to", outcome_trait),
                  exposure_pop, outcome_pop,
                  paste0(round(mismatch_estimate, 3), " (", mismatch_p, ")"),
                  paste0(round(reference_estimate, 3), " (", reference_p, ")"),
                  zconcord),
       col.names = FALSE,
       "reproducible_figures/tables4_grapple_mismatchlarger.csv")

#=Supplemental Table S5 - same results applied to IVW and BEE=
mtr.ratio.bee1e5 <- vector()
reference.estimate.bee1e5 <- vector()
reference.se.bee1e5 <- vector()

mtr.ratio.bee5e8 <- vector()
reference.estimate.bee5e8 <- vector()
reference.se.bee5e8 <- vector()

mtr.ratio.ivw <- vector()
reference.estimate.ivw <- vector()
reference.se.ivw <- vector()

for (i in 1:dim(mr_table.exactETOS)[1]) {
  if (mr_table.exactETOS$is.bestmatch[i]) { #mismatch/reference ratio only applies to mismatches.
    mtr.ratio.bee1e5[i] <- NA
    reference.estimate.bee1e5[i] <- NA
    reference.se.bee1e5[i] <- NA
    mtr.ratio.bee5e8[i] <- NA
    reference.estimate.bee5e8[i] <- NA
    reference.se.bee5e8[i] <- NA
    mtr.ratio.ivw[i] <- NA
    reference.estimate.ivw[i] <- NA
    reference.se.ivw[i] <- NA
    next
  }

  mismatch.estimate.bee1e5 <- mr_table.exactETOS$bee1e5.beta.rescale[i]
  mismatch.estimate.bee5e8 <- mr_table.exactETOS$bee5e8.beta.rescale[i]
  mismatch.estimate.ivw <- mr_table.exactETOS$ivw.beta.rescale[i]
  ETOS <- mr_table.exactETOS$ETOS[i]

  reference.estimate.bee1e5[i] <- mr_table.exactETOS$bee1e5.beta.rescale[mr_table.exactETOS$is.bestmatch & mr_table.exactETOS$ETOS == ETOS]
  reference.se.bee1e5[i] <- mr_table.exactETOS$bee1e5.se.rescale[mr_table.exactETOS$is.bestmatch & mr_table.exactETOS$ETOS == ETOS]
  mtr.ratio.bee1e5[i] <- mismatch.estimate.bee1e5/reference.estimate.bee1e5[i]

  if (!is.na(mismatch.estimate.bee5e8)) {
    reference.estimate.bee5e8[i] <- mr_table.exactETOS$bee5e8.beta.rescale[mr_table.exactETOS$is.bestmatch & mr_table.exactETOS$ETOS == ETOS]
    reference.se.bee5e8[i] <- mr_table.exactETOS$bee5e8.se.rescale[mr_table.exactETOS$is.bestmatch & mr_table.exactETOS$ETOS == ETOS]
    mtr.ratio.bee5e8[i] <- mismatch.estimate.bee5e8/reference.estimate.bee5e8[i]
    reference.estimate.ivw[i] <- mr_table.exactETOS$ivw.beta.rescale[mr_table.exactETOS$is.bestmatch & mr_table.exactETOS$ETOS == ETOS]
    reference.se.ivw[i] <- mr_table.exactETOS$ivw.se.rescale[mr_table.exactETOS$is.bestmatch & mr_table.exactETOS$ETOS == ETOS]
    mtr.ratio.ivw[i] <- mismatch.estimate.ivw/reference.estimate.ivw[i]
  } else {
    reference.estimate.bee5e8[i] <- NA
    reference.se.bee5e8[i] <- NA
    mtr.ratio.bee5e8[i] <- NA
    reference.estimate.ivw[i] <- NA
    reference.se.ivw[i] <- NA
    mtr.ratio.ivw[i] <- NA
  }
}

#Calculating the overall shrinkage coefficients for IVW/BEE using SIMEX
set.seed(2025)
indices.1e5 <- which(!is.na(mr_table.exactETOS$zconcord.bee1e5))
simex.bee1e5 <- simex_regression(outcome = mr_table.exactETOS$bee1e5.beta.rescale[indices.1e5],
                                     outcome.se = mr_table.exactETOS$bee1e5.se.rescale[indices.1e5],
                                     predictor = reference.estimate.1e5[indices.1e5],
                                     predictor.se = reference.se.1e5[indices.1e5])

indices.5e8 <- which(!is.na(mr_table.exactETOS$zconcord.bee5e8))
simex.bee5e8 <- simex_regression(outcome = mr_table.exactETOS$bee5e8.beta.rescale[indices.5e8],
                                     outcome.se = mr_table.exactETOS$bee5e8.se.rescale[indices.5e8],
                                     predictor = reference.estimate.5e8[indices.5e8],
                                     predictor.se = reference.se.5e8[indices.5e8])
simex.ivw <- simex_regression(outcome = mr_table.exactETOS$ivw.beta.rescale[indices.5e8],
                                     outcome.se = mr_table.exactETOS$ivw.se.rescale[indices.5e8],
                                     predictor = reference.estimate.5e8[indices.5e8],
                                     predictor.se = reference.se.5e8[indices.5e8])

table_categories <- c("MR Method",
                      "Total Number of Mismatching Estimates with Reference Estimates",
                      "Mismatching Estimates with Non-Significant Concordance z-scores",
                      "Mismatching Estimates with Significant Concordance z-scores, Same Sign as Reference, Smaller Magnitude",
                      "Mismatching Estimates with Significant Concordance z-scores, Same Sign as Reference, Larger Magnitude",
                      "Mismatching Estimates with Significant Concordance z-scores, Opposite Sign as Reference",
                      "Mismatching Estimates with Significant Concordance z-scores, Opposite Sign as Reference, Reference and Mismatch MR p-values < 0.05",
                      "SIMEX-Estimated Overall Shrinkage Coefficient, 1000 Bootstrap Resamples")

numbers_ivw <- c("IVW",
                 sum(!is.na(mr_table.exactETOS$zconcord.ivw)),
                 sum(mr_table.exactETOS$pconcord.ivw.bh >= 0.05, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.ivw.bh < 0.05 & mtr.ratio.5e8 > 0 & mtr.ratio.5e8 < 1, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.ivw.bh < 0.05 & mtr.ratio.5e8 > 1, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.ivw.bh < 0.05 & mtr.ratio.5e8 < 0, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.ivw.bh < 0.05 & mtr.ratio.5e8 < 0 &
                     abs(mr_table.exactETOS$ivw.beta.rescale/mr_table.exactETOS$ivw.se.rescale) > crit.z &
                     abs(reference.estimate.ivw/reference.se.ivw) > crit.z, na.rm = TRUE),
                 paste0(round(mean(simex.ivw), 3), " (", round(sd(simex.ivw), 3), ")"))

numbers_bee1e5 <- c("BEE (1e-5)",
                 sum(!is.na(mr_table.exactETOS$zconcord.bee1e5)),
                 sum(mr_table.exactETOS$pconcord.bee1e5.bh >= 0.05, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.bee1e5.bh < 0.05 & mtr.ratio.bee1e5 > 0 & mtr.ratio.bee1e5 < 1, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.bee1e5.bh < 0.05 & mtr.ratio.bee1e5 > 1, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.bee1e5.bh < 0.05 & mtr.ratio.bee1e5 < 0, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.bee1e5.bh < 0.05 & mtr.ratio.bee1e5 < 0 &
                     abs(mr_table.exactETOS$bee1e5.beta.rescale/mr_table.exactETOS$bee1e5.se.rescale) > crit.z &
                     abs(reference.estimate.bee1e5/reference.se.bee1e5) > crit.z, na.rm = TRUE),
                 paste0(round(mean(simex.bee1e5), 3), " (", round(sd(simex.bee1e5), 3), ")"))

numbers_bee5e8 <- c("BEE (5e-8)",
                 sum(!is.na(mr_table.exactETOS$zconcord.bee5e8)),
                 sum(mr_table.exactETOS$pconcord.bee5e8.bh >= 0.05, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.bee5e8.bh < 0.05 & mtr.ratio.5e8 > 0 & mtr.ratio.5e8 < 1, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.bee5e8.bh < 0.05 & mtr.ratio.5e8 > 1, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.bee5e8.bh < 0.05 & mtr.ratio.5e8 < 0, na.rm = TRUE),
                 sum(mr_table.exactETOS$pconcord.bee5e8.bh < 0.05 & mtr.ratio.5e8 < 0 &
                     abs(mr_table.exactETOS$bee5e8.beta.rescale/mr_table.exactETOS$bee5e8.se.rescale) > crit.z &
                     abs(reference.estimate.bee5e8/reference.se.bee5e8) > crit.z, na.rm = TRUE),
                 paste0(round(mean(simex.bee5e8), 3), " (", round(sd(simex.bee5e8), 3), ")"))

fwrite(data.frame(table_categories, numbers_ivw, numbers_bee1e5, numbers_5e8), col.names = FALSE,
       "reproducible_figures/tables5_bee_ivw.csv")

#==Subpopulation-specific SIMEX + Fst analyses==
source("simex_source_func.R")

#Incorporating population information
grapple1e5_pops <- fread("study_info/pop_id_reference.csv")
mr_table.exactETOS.pops <- left_join(mr_table.exactETOS, grapple1e5_pops, by = c("exposure_id", "outcome_id"))

subpopulations <- intersect(mr_table.exactETOS.pops$exposure_pop, mr_table.exactETOS.pops$outcome_pop)
mr_list.ETOS <- split(mr_table.exactETOS.pops, f = mr_table.exactETOS.pops$ETOS)

set.seed(2025)
master_list <- list()
counter <- 0

for (i in 1:length(subpopulations)) {
  target_pop <- subpopulations[i]
  selection_pops <- subpopulations[-i]
  for (selection_pop in selection_pops) {
    counter <- counter + 1
    print(paste("=====Starting", selection_pop, "to", paste0(target_pop,"=====")))
    ETOS.pair.selection <- lapply(mr_list.ETOS, function(x){return(unique(x$outcome_pop) == target_pop &
                                                                     any(x$exposure_pop == target_pop) &
                                                                     any(x$exposure_pop == selection_pop))}) %>% unlist()
    ETOS.pairs <- mr_list.ETOS[ETOS.pair.selection]

    ETOS.pair <- lapply(ETOS.pairs, function(x){return(unique(x$ETOS))}) %>% unlist() %>% unname()
    ETOS.pair.target_outcome <- lapply(ETOS.pairs, function(x){return(unique(x$outcome_id))}) %>% unlist() %>% unname()

    ETOS.pair.target_exposure <- lapply(ETOS.pairs, function(x){targets <- x[x$exposure_pop == target_pop,]; targets2 <- targets[which.max(targets$exposure.n),]; return(targets2$exposure_id)}) %>% unlist() %>% unname()
    ETOS.pair.selection_exposure <- lapply(ETOS.pairs, function(x){targets <- x[x$exposure_pop == selection_pop,]; targets2 <- targets[which.max(targets$exposure.n),]; return(targets2$exposure_id)}) %>% unlist %>% unname()

    #==Estimates for GRAPPLE, 1e-5==
    target_to_target.beta_1e5 <- lapply(ETOS.pairs, function(x){targets <- x[x$exposure_pop == target_pop,]; targets2 <- targets[which.max(targets$exposure.n),]; return(targets2$grapple1e5.beta.rescale)}) %>% unlist() %>% unname()
    target_to_target.se_1e5 <- lapply(ETOS.pairs, function(x){targets <- x[x$exposure_pop == target_pop,]; targets2 <- targets[which.max(targets$exposure.n),]; return(targets2$grapple1e5.se.rescale)}) %>% unlist() %>% unname()

    selection_to_target.beta_1e5 <- lapply(ETOS.pairs, function(x){targets <- x[x$exposure_pop == selection_pop,]; targets2 <- targets[which.max(targets$exposure.n),]; return(targets2$grapple1e5.beta.rescale)}) %>% unlist %>% unname()
    selection_to_target.se_1e5 <- lapply(ETOS.pairs, function(x){targets <- x[x$exposure_pop == selection_pop,]; targets2 <- targets[which.max(targets$exposure.n),]; return(targets2$grapple1e5.se.rescale)}) %>% unlist %>% unname()

    #==Estimates for GRAPPLE, 5e-8==
    target_to_target.beta_5e8 <- lapply(ETOS.pairs, function(x){targets <- x[x$exposure_pop == target_pop,]; targets2 <- targets[which.max(targets$exposure.n),]; return(targets2$grapple5e8.beta.rescale)}) %>% unlist() %>% unname()
    target_to_target.se_5e8 <- lapply(ETOS.pairs, function(x){targets <- x[x$exposure_pop == target_pop,]; targets2 <- targets[which.max(targets$exposure.n),]; return(targets2$grapple5e8.se.rescale)}) %>% unlist() %>% unname()

    selection_to_target.beta_5e8 <- lapply(ETOS.pairs, function(x){targets <- x[x$exposure_pop == selection_pop,]; targets2 <- targets[which.max(targets$exposure.n),]; return(targets2$grapple5e8.beta.rescale)}) %>% unlist %>% unname()
    selection_to_target.se_5e8 <- lapply(ETOS.pairs, function(x){targets <- x[x$exposure_pop == selection_pop,]; targets2 <- targets[which.max(targets$exposure.n),]; return(targets2$grapple5e8.se.rescale)}) %>% unlist %>% unname()

    ETOS.table <- data.frame(ETOS.pair, outcome.id = ETOS.pair.target_outcome,
                             target_exposure.id = ETOS.pair.target_exposure, selection_exposure.id = ETOS.pair.selection_exposure,
                             samepop.1e5 = target_to_target.beta_1e5, samepop.se.1e5 = target_to_target.se_1e5,
                             crosspop.1e5 = selection_to_target.beta_1e5, crosspop.se.1e5 = selection_to_target.se_1e5,
                             samepop.5e8 = target_to_target.beta_5e8, samepop.se.5e8 = target_to_target.se_5e8,
                             crosspop.5e8 = selection_to_target.beta_5e8, crosspop.se.5e8 = selection_to_target.se_5e8)

    if (dim(ETOS.table)[1] > 0) {
      ETOS.table.1e5 <- ETOS.table[,c(1:8)]
      ETOS.table.1e5 <- ETOS.table.1e5[complete.cases(ETOS.table.1e5),]

      ETOS.table.5e8 <- ETOS.table[,c(1:4, 9:12)]
      ETOS.table.5e8 <- ETOS.table.5e8[complete.cases(ETOS.table.5e8),]
    } else { #zero-dimensional case
      ETOS.table.1e5 <- ETOS.table
      ETOS.table.5e8 <- ETOS.table
    }

    #Now, we obtain our averages via SIMEX regression.
    print("==1e-5==")
    if (dim(ETOS.table.1e5)[1] >= 8) {
      simex_resamples.1e5 <- simex_regression(predictor = ETOS.table.1e5$samepop.1e5, predictor.se = ETOS.table.1e5$samepop.se.1e5,
                         outcome = ETOS.table.1e5$crosspop.1e5, outcome.se = ETOS.table.1e5$crosspop.se.1e5)
      simex_mean.1e5 <- mean(simex_resamples.1e5)
      simex_se.1e5 <- sd(simex_resamples.1e5)
    } else {
      print("Not at least 8 ET-OS pairs, skipping SIMEX regression for 1e-5.")
      simex_resamples.1e5 <- "Insufficient ET-OS pairs"
      simex_mean.1e5 <- "Insufficient ET-OS pairs"
      simex_se.1e5 <- "Insufficient ET-OS pairs"
    }

    print("==5e-8==")
    if (dim(ETOS.table.5e8)[1] >= 8) {
      simex_resamples.5e8 <- simex_regression(predictor = ETOS.table.5e8$samepop.5e8, predictor.se = ETOS.table.5e8$samepop.se.5e8,
                         outcome = ETOS.table.5e8$crosspop.5e8, outcome.se = ETOS.table.5e8$crosspop.se.5e8)
      simex_mean.5e8 <- mean(simex_resamples.5e8)
      simex_se.5e8 <- sd(simex_resamples.5e8)
    } else {
      print("Not at least 8 ET-OS pairs, skipping SIMEX regression for 5e-8.")
      simex_resamples.5e8 <- "Insufficient ET-OS pairs"
      simex_mean.5e8 <- "Insufficient ET-OS pairs"
      simex_se.5e8 <- "Insufficient ET-OS pairs"
    }

    master_list[[counter]] <- list(estimate_table = ETOS.table,
                                   simex_resamples.1e5 = simex_resamples.1e5, simex_mean.1e5 = simex_mean.1e5, simex_se.1e5 = simex_se.1e5,
                                   simex_resamples.5e8 = simex_resamples.5e8, simex_mean.5e8 = simex_mean.5e8, simex_se.5e8 = simex_se.5e8)
    names(master_list)[counter] <- paste(selection_pop, "to", target_pop)
  }
}

#saveRDS(master_list, "/net/snowwhite/home/lijack/6_10_2024_MR_survey/2025_12_18/master_list_simex.RDS")

pop_pair <- names(master_list)
simex.1e5 <- master_list %>% map("simex_mean.1e5") %>% unlist()
simex_se.1e5 <- master_list %>% map("simex_se.1e5") %>% unlist()
simex.5e8 <- master_list %>% map("simex_mean.5e8") %>% unlist()
simex_se.5e8 <- master_list %>% map("simex_se.5e8") %>% unlist()
ETOS.pair.count.1e5 <- lapply(master_list, function(x){return(dim(x$estimate_table)[1])}) %>% unlist()
ETOS.pair.count.5e8 <- lapply(master_list, function(x){if (dim(x$estimate_table)[1] == 0) {return(0)} else {ETOS.table.5e8 <- x$estimate_table[,c(1:4, 9:12)]; ETOS.table.5e8 <- ETOS.table.5e8[complete.cases(ETOS.table.5e8),]; return(dim(ETOS.table.5e8)[1])}}) %>% unlist()

#subpopulation order
#[1] "Biobank Japan"
#[2] "SAS (Pan-UK Biobank)"
#[3] "UK Biobank"
#[4] "AFR (Pan-UK Biobank)"
#[5] "Majority Hispanic mega-analysis (Wojcik)"
#[6] "FinnGen (Release 10)"
#[7] "Taiwan Biobank"

#=Supplemental Figure S2 - Matrices of SIMEX slopes=
selection <- gsub("(.*) to (.*)", "\\1", pop_pair)
target <- gsub("(.*) to (.*)", "\\2", pop_pair)

simex1e5.matrix <- matrix(NA, nrow = 7, ncol = 7)
simex5e8.matrix <- matrix(NA, nrow = 7, ncol = 7)

rownames(simex1e5.matrix) <- colnames(simex1e5.matrix) <- subpopulations[c(3, 6, 1, 7, 5, 2, 4)]
rownames(simex5e8.matrix) <- colnames(simex5e8.matrix) <- subpopulations[c(3, 6, 1, 7, 5, 2, 4)]

for (i in 1:length(master_list)) {
  row_index <- which(rownames(simex1e5.matrix) == selection[i])
  col_index <- which(colnames(simex1e5.matrix) == target[i])

  if (simex.1e5[i] == "Insufficient ET-OS pairs") {
    simex1e5.matrix[row_index, col_index] <- NA
  } else {
    simex1e5.matrix[row_index, col_index] <- paste0(round(as.numeric(simex.1e5[i]), 3),
                                                    " (", round(as.numeric(simex_se.1e5[i]), 3),
                                                    ") [", ETOS.pair.count.1e5[i], "]")
  }

  if (simex.5e8[i] == "Insufficient ET-OS pairs") {
    simex5e8.matrix[row_index, col_index] <- NA
  } else {
    simex5e8.matrix[row_index, col_index] <- paste0(round(as.numeric(simex.5e8[i]), 3),
                                                    " (", round(as.numeric(simex_se.5e8[i]), 3),
                                                    ") [", ETOS.pair.count.5e8[i], "]")
  }
}

write.csv(simex1e5.matrix, "reproducible_figures/figures2_simex1e5.csv")
write.csv(simex5e8.matrix, "reproducible_figures/figures2_simex5e8.csv")

fixation_indices <- c(0.083, 0.118, 0.212, NA, 0.110, 0.009,
                      0.083, 0.019, 0.156, NA, 0.021, 0.085,
                      0.118, 0.019, 0.178, NA, 0.007, 0.119,
                      0.212, 0.156, 0.178, NA, 0.179, 0.213,
                      NA, NA, NA, NA, NA, NA,
                      0.110, 0.021, 0.007, 0.179, NA, 0.111,
                      0.009, 0.085, 0.119, 0.213, NA, 0.111)

simex_table.1e5 <- data.frame(pop_pair = pop_pair,
                              simex = simex.1e5 %>% as.numeric(), simex_se = simex_se.1e5 %>% as.numeric(),
                              ETOS.pair.count = ETOS.pair.count.1e5, fst = fixation_indices)
simex_table.1e5 <- simex_table.1e5[complete.cases(simex_table.1e5[,1:4]),]

simex_table.5e8 <- data.frame(pop_pair = pop_pair,
                              simex = simex.5e8 %>% as.numeric(), simex_se = simex_se.5e8 %>% as.numeric(),
                              ETOS.pair.count = ETOS.pair.count.5e8, fst = fixation_indices)
simex_table.5e8 <- simex_table.5e8[complete.cases(simex_table.5e8[,1:4]),]

#=Figure 4 - Fst scatterplot=

summary(lm(simex_table.1e5$simex ~ simex_table.1e5$fst, weights = (1/simex_table.1e5$simex_se)^2))
fst.slope.1e5 <- lm(simex_table.1e5$simex ~ simex_table.1e5$fst, weights = (1/simex_table.1e5$simex_se)^2)

summary(lm(simex_table.5e8$simex ~ simex_table.5e8$fst, weights = (1/simex_table.5e8$simex_se)^2))
fst.slope.5e8 <- lm(simex_table.5e8$simex ~ simex_table.5e8$fst, weights = (1/simex_table.5e8$simex_se)^2)

ggplot(data = simex_table.1e5, aes(x = fst, y = simex, col = factor(fst))) + geom_point(size = 2) + labs(title = "", x = "Approximate Fixation Index (Fst)", y = "Estimated Cross- vs. Same-Population\nRegression Slope") + theme_bw(base_size = 17.5) + geom_hline(yintercept = 0) + geom_abline(slope = fst.slope.1e5$coefficients[2], intercept = fst.slope.1e5$coefficients[1], lty = 2) + geom_errorbar(mapping = aes(ymin = simex - 1.96*simex_se, ymax = simex + 1.96*simex_se)) + guides(col ="none")
ggsave("reproducible_figures/figure4_fst1e5.png")

ggplot(data = simex_table.5e8, aes(x = fst, y = simex, col = factor(fst))) + geom_point(size = 2) + labs(title = "", x = "Approximate Fixation Index (Fst)", y = "Estimated Cross- vs. Same-Population\nRegression Slope") + theme_bw(base_size = 17.5) + geom_hline(yintercept = 0) + geom_abline(slope = fst.slope.5e8$coefficients[2], intercept = fst.slope.5e8$coefficients[1], lty = 2) + geom_errorbar(mapping = aes(ymin = simex - 1.96*simex_se, ymax = simex + 1.96*simex_se)) + guides(col ="none")
ggsave("reproducible_figures/figure4_fst5e8.png")

#=Supplemental Figure S4 - Scatterplots showing SIMEX slopes for specific subpopulations=

#shortening subpopulation names
shortened.names <- function(x) {
  if (x == "Majority Hispanic mega-analysis (Wojcik)") {
    return("Majority Hispanic")
  }
  if (x == "SAS (Pan-UK Biobank)") {
    return("SAS (Pan-UKB)")
  }
  if (x == "FinnGen (Release 10)") {
    return("FinnGen")
  }
  if (x == "AFR (Pan-UK Biobank)") {
    return("AFR (Pan-UKB)")
  } else {
    return(x)
  }
}


for (pair in intersect(simex_table.1e5$pop_pair, simex_table.5e8$pop_pair)) {
  i <- which(names(master_list) == pair)
  target.pop <- gsub("(.*) to (.*)", "\\2", pair) %>% shortened.names()
  selection.pop <- gsub("(.*) to (.*)", "\\1", pair) %>% shortened.names()

  png(file = paste0("reproducible_figures/figure_s4/", target.pop, "_", selection.pop,".png"),
      width = 1920, height = 1080)
  par(mar = c(5.1, 5.1, 4.1, 2.1), mfcol = c(1,2))

  plot.min <- min(min(master_list[[i]]$estimate_table$samepop.1e5, na.rm = TRUE) - 0.02,
                  min(master_list[[i]]$estimate_table$crosspop.1e5, na.rm = TRUE) - 0.02, na.rm = TRUE)
  plot.max <- max(max(master_list[[i]]$estimate_table$samepop.1e5, na.rm = TRUE) + 0.02,
                  max(master_list[[i]]$estimate_table$crosspop.1e5, na.rm = TRUE) + 0.02, na.rm = TRUE)

  plot(master_list[[i]]$estimate_table$samepop.1e5, master_list[[i]]$estimate_table$crosspop.1e5,
     xlim = c(plot.min, plot.max),
     ylim = c(plot.min, plot.max),
     xlab = paste0(target.pop, " to ", target.pop),
     ylab = paste0(selection.pop, " to ", target.pop),
     main = "GRAPPLE (1e-5)",
     cex.main = 2, cex.lab = 2, cex.axis = 1.5, pch = 20)
  abline(h = 0, col = "black")
  abline(v = 0, col = "black")
  abline(a = 0, b = 1, col = "black")
  abline(a = 0, b = master_list[[i]]$simex_mean.1e5, lty = 2)

  plot.min <- min(min(master_list[[i]]$estimate_table$samepop.5e8) - 0.02,
                  min(master_list[[i]]$estimate_table$crosspop.5e8, na.rm = TRUE) - 0.02, na.rm = TRUE)
  plot.max <- max(max(master_list[[i]]$estimate_table$samepop.5e8) + 0.02,
                  max(master_list[[i]]$estimate_table$crosspop.5e8, na.rm = TRUE) + 0.02, na.rm = TRUE)


  plot(master_list[[i]]$estimate_table$samepop.5e8, master_list[[i]]$estimate_table$crosspop.5e8,
     xlim = c(plot.min, plot.max),
     ylim = c(plot.min, plot.max),
     xlab = paste0(target.pop, " to ", target.pop),
     ylab = paste0(selection.pop, " to ", target.pop),
     main = "GRAPPLE (5e-8)",
     cex.main = 2, cex.lab = 2, cex.axis = 1.5, pch = 20)
  abline(h = 0, col = "black")
  abline(v = 0, col = "black")
  abline(a = 0, b = 1, col = "black")
  abline(a = 0, b = master_list[[i]]$simex_mean.5e8, lty = 2)

  dev.off()
  print(pair)
}

#that leaves just the population pairs specific to 1e-5, most of which are African
  #We can just use the same structure, I think

for (pair in setdiff(simex_table.1e5$pop_pair, simex_table.5e8$pop_pair)) {
  i <- which(names(master_list) == pair)
  target.pop <- gsub("(.*) to (.*)", "\\2", pair) %>% shortened.names()
  selection.pop <- gsub("(.*) to (.*)", "\\1", pair) %>% shortened.names()

  png(file = paste0("reproducible_figures/figure_s4/", target.pop, "_", selection.pop,".png"),
      width = 1920, height = 1080)
  par(mar = c(5.1, 5.1, 4.1, 2.1), mfcol = c(1,2))
  plot.min <- min(min(master_list[[i]]$estimate_table$samepop.1e5) - 0.02, min(master_list[[i]]$estimate_table$crosspop.1e5) - 0.02)
  plot.max <- max(max(master_list[[i]]$estimate_table$samepop.1e5) + 0.02, max(master_list[[i]]$estimate_table$crosspop.1e5) + 0.02)

  plot(master_list[[i]]$estimate_table$samepop.1e5, master_list[[i]]$estimate_table$crosspop.1e5,
     xlim = c(plot.min, plot.max),
     ylim = c(plot.min, plot.max),
     xlab = paste0(target.pop, " to ", target.pop),
     ylab = paste0(selection.pop, " to ", target.pop),
     main = "GRAPPLE (1e-5)",
     cex.main = 2, cex.lab = 2, cex.axis = 1.5, pch = 20)
  abline(h = 0, col = "black")
  abline(v = 0, col = "black")
  abline(a = 0, b = 1, col = "black")
  abline(a = 0, b = master_list[[i]]$simex_mean.1e5, lty = 2)

  dev.off()
  print(pair)
}

#==UKBB/BBJ to BBJ and replication analysis==
mr_table.pops <- left_join(mr_table, grapple1e5_pops, by = c("exposure_id", "outcome_id"))

exposure_info <- fread("study_info/exposure_info.csv")
outcome_info <- fread("study_info/outcome_info.csv")

mr_table.pops$exposure_trait <- sapply(mr_table.pops$exposure_id, function(X){exposure_info$Trait[exposure_info$`Study ID` == X]}) %>% unlist() %>% unname()
mr_table.pops$outcome_trait <- sapply(mr_table.pops$outcome_id, function(X){outcome_info$Trait[outcome_info$`Study ID` == X]}) %>% unlist() %>% unname()

mr_table.pops$trait_pair <- paste0(mr_table.pops$exposure_trait, "-", mr_table.pops$outcome_trait)

bbj_to_bbj_traitpairs <- unique(mr_table.pops$trait_pair[mr_table.pops$exposure_pop == "Biobank Japan" & mr_table.pops$outcome_pop == "Biobank Japan"])
sakaue_to_bbj_traitpairs <- unique(mr_table.pops$trait_pair[mr_table.pops$exposure_pop == "UKBB + BBJ meta-analysis" & mr_table.pops$outcome_pop == "Biobank Japan"])

traitpairs <- intersect(bbj_to_bbj_traitpairs, sakaue_to_bbj_traitpairs)

outcome.id <- vector()
samepop.id <- vector()
crosspop.id <- vector()

samepop.1e5 <- vector()
samepop.se.1e5 <- vector()
samepop.p.1e5 <- vector()
crosspop.1e5 <- vector()
crosspop.se.1e5 <- vector()
crosspop.p.1e5 <- vector()

samepop.5e8 <- vector()
samepop.se.5e8 <- vector()
samepop.p.5e8 <- vector()
crosspop.5e8 <- vector()
crosspop.se.5e8 <- vector()
crosspop.p.5e8 <- vector()

for (i in 1:length(traitpairs)) {
  traitpair_estimates <- mr_table.pops[mr_table.pops$trait_pair == traitpairs[i],]
  samepop_index <- which(traitpair_estimates$exposure_pop == "Biobank Japan" & traitpair_estimates$outcome_pop == "Biobank Japan")
  crosspop_index <- which(traitpair_estimates$exposure_pop == "UKBB + BBJ meta-analysis" & traitpair_estimates$outcome_pop == "Biobank Japan")

  outcome.id[i] <- traitpair_estimates$outcome_id[samepop_index]
  samepop.id[i] <- traitpair_estimates$exposure_id[samepop_index]
  crosspop.id[i] <- traitpair_estimates$exposure_id[crosspop_index]

  samepop.1e5[i] <- traitpair_estimates$grapple1e5.beta[samepop_index]
  samepop.se.1e5[i] <- traitpair_estimates$grapple1e5.se[samepop_index]
  samepop.p.1e5[i] <- 2 * pnorm(abs(samepop.1e5[i]/samepop.se.1e5[i]), lower.tail = FALSE)
  crosspop.1e5[i] <- traitpair_estimates$grapple1e5.beta[crosspop_index]
  crosspop.se.1e5[i] <- traitpair_estimates$grapple1e5.se[crosspop_index]
  crosspop.p.1e5[i] <- 2 * pnorm(abs(crosspop.1e5[i]/crosspop.se.1e5[i]), lower.tail = FALSE)

  samepop.5e8[i] <- traitpair_estimates$grapple5e8.beta[samepop_index]
  samepop.se.5e8[i] <- traitpair_estimates$grapple5e8.se[samepop_index]
  samepop.p.5e8[i] <- 2 * pnorm(abs(samepop.5e8[i]/samepop.se.5e8[i]), lower.tail = FALSE)
  crosspop.5e8[i] <- traitpair_estimates$grapple5e8.beta[crosspop_index]
  crosspop.se.5e8[i] <- traitpair_estimates$grapple5e8.se[crosspop_index]
  crosspop.p.5e8[i] <- 2 * pnorm(abs(crosspop.5e8[i]/crosspop.se.5e8[i]), lower.tail = FALSE)

}

sakaue_bbj_table <- data.frame(traitpairs, outcome.id, samepop.id, crosspop.id,
                                      samepop.1e5, samepop.se.1e5, samepop.p.1e5,
                                      crosspop.1e5, crosspop.se.1e5, crosspop.p.1e5,
                                      samepop.5e8, samepop.se.5e8, samepop.p.5e8,
                                      crosspop.5e8, crosspop.se.5e8, crosspop.p.5e8)

sakaue_bbj_table$sig_status.1e5 <- ifelse(sakaue_bbj_table$crosspop.p.1e5 < 0.05/195,
                                      ifelse(sakaue_bbj_table$samepop.p.1e5 < 0.05/195, "Significance Maintained", "Gain of Significance"),
                                      ifelse(sakaue_bbj_table$samepop.p.1e5 < 0.05/195, "Loss of Significance", "Non-Significance Maintained"))

sakaue_bbj_table$sig_status.5e8 <- ifelse(sakaue_bbj_table$crosspop.p.5e8 < 0.05/195,
                                      ifelse(sakaue_bbj_table$samepop.p.5e8 < 0.05/195, "Significance Maintained", "Gain of Significance"),
                                      ifelse(sakaue_bbj_table$samepop.p.5e8 < 0.05/195, "Loss of Significance", "Non-Significance Maintained"))

#the trait pairs significant for ONLY ONE of same-pop or cross-pop at a threshold of 1e-5
#sig_ukbb_bbj_1e5 <- which(sakaue_bbj_table$crosspop.p.1e5 < 0.05/195 | sakaue_bbj_table$samepop.p.1e5 < 0.05/195)
sig_ukbb_bbj_1e5 <- which((sakaue_bbj_table$crosspop.p.1e5 < 0.05/195 & sakaue_bbj_table$samepop.p.1e5 > 0.05/195) |
                          (sakaue_bbj_table$crosspop.p.1e5 > 0.05/195 & sakaue_bbj_table$samepop.p.1e5 < 0.05/195))
length(sig_ukbb_bbj_1e5)

#sig_ukbb_bbj_5e8 <- which(sakaue_bbj_table$crosspop.p.5e8 < 0.05/195 | sakaue_bbj_table$samepop.p.5e8 < 0.05/195)
sig_ukbb_bbj_5e8 <- which((sakaue_bbj_table$crosspop.p.5e8 < 0.05/195 & sakaue_bbj_table$samepop.p.5e8 > 0.05/195) |
                          (sakaue_bbj_table$crosspop.p.5e8 > 0.05/195 & sakaue_bbj_table$samepop.p.5e8 < 0.05/195))
length(sig_ukbb_bbj_5e8)

#Preparing summary tables of MR estimates for replication

#Replication at 1e-5
replication_table.1e5 <- sakaue_bbj_table[sig_ukbb_bbj_1e5,] %>% select(!which(names(sakaue_bbj_table) %like% "5e8"))
rep.pairs.available <- vector() #If that trait pair has any non-UK, non-Japan exact-matching 2SMR estimates at a threshold of 1e-5

rep1.beta <- vector()
rep1.pvalue <- vector()
rep1.targetpop <- vector()

rep2.beta <- vector()
rep2.pvalue <- vector()
rep2.targetpop <- vector()

rep.minp <- vector()

for (i in 1:length(sig_ukbb_bbj_1e5)) {
  traitpair.estimates.1e5 <- mr_table.pops[mr_table.pops$trait_pair == replication_table.1e5$traitpairs[i] &
                                                   !(mr_table.pops$exposure_pop %in% c("UK Biobank", "Biobank Japan")) &
                                                   mr_table.pops$exposure_pop == mr_table.pops$outcome_pop,] %>%
                                                   select("exposure_pop", "grapple1e5.beta", "grapple1e5.se")
  traitpair.estimates.1e5$grapple1e5.p <- 2 * pnorm(abs(traitpair.estimates.1e5$grapple1e5.beta/traitpair.estimates.1e5$grapple1e5.se),
                                                    lower.tail = FALSE)

  if (dim(traitpair.estimates.1e5)[1] >=1 ) {
    rep.pairs.available[i] <- dim(traitpair.estimates.1e5)[1]

    #we grab all summary statistics
    rep1.beta[i] <- traitpair.estimates.1e5$grapple1e5.beta[1]
    rep1.pvalue[i] <- traitpair.estimates.1e5$grapple1e5.p[1]
    rep1.targetpop[i] <- traitpair.estimates.1e5$exposure_pop[1]

    rep2.beta[i] <- NA
    rep2.pvalue[i] <- NA
    rep2.targetpop[i] <- NA

    #we have only 1 or 2 pairs available for replication, which means I can get away with some pretty crappy code
    if (rep.pairs.available[i] == 2) {
        rep2.beta[i] <- traitpair.estimates.1e5$grapple1e5.beta[2]
        rep2.pvalue[i] <- traitpair.estimates.1e5$grapple1e5.p[2]
        rep2.targetpop[i] <- traitpair.estimates.1e5$exposure_pop[2]
    }

    sign_align.1e5 <- traitpair.estimates.1e5[sign(traitpair.estimates.1e5$grapple1e5.beta) == sign(replication_table.1e5$crosspop.1e5[i]),]
    if (dim(sign_align.1e5)[1] >= 1) {
      rep.minp[i] <- min(sign_align.1e5$grapple1e5.p)
    } else { #there are no estimates with the same sign as the cross-population estimate.
      rep.minp[i] <- 1
    }
  } else {
    rep.pairs.available[i] <- 0
    rep1.beta[i] <- rep2.beta[i] <- rep1.pvalue[i] <- rep2.pvalue[i] <- rep1.targetpop[i] <- rep2.targetpop[i] <- NA
    rep.minp[i] <- NA
  }
}

replication_table.1e5$rep.pairs.available <- rep.pairs.available
replication_table.1e5$rep1.beta <- rep1.beta
replication_table.1e5$rep1.pvalue <- rep1.pvalue
replication_table.1e5$rep1.targetpop <- rep1.targetpop
replication_table.1e5$rep2.beta <- rep2.beta
replication_table.1e5$rep2.pvalue <- rep2.pvalue
replication_table.1e5$rep.minp <- rep.minp
replication_table.1e5$rep2.targetpop <- rep2.targetpop

rep.threshold.1e5 <- 0.05/sum(replication_table.1e5$rep.pairs.available > 0)
replication_table.1e5$replicated.1e5 <- ifelse(replication_table.1e5$rep.pairs.available > 0,
                                           ifelse(replication_table.1e5$rep.minp < rep.threshold.1e5, "Replicated", "Not Replicated"),
                                           "Not Tested")


#Replication at 5e-8
replication_table.5e8 <- sakaue_bbj_table[sig_ukbb_bbj_5e8,] %>% select(!which(names(sakaue_bbj_table) %like% "1e5"))
rep.pairs.available <- vector() #If that trait pair has any non-UK, non-Japan exact-matching 2SMR estimates at a threshold of 1e-5

rep1.beta <- vector()
rep1.pvalue <- vector()
rep1.targetpop <- vector()

rep2.beta <- vector()
rep2.pvalue <- vector()
rep2.targetpop <- vector()

rep.minp <- vector()

for (i in 1:length(sig_ukbb_bbj_5e8)) {
  traitpair.estimates.5e8 <- mr_table.pops[mr_table.pops$trait_pair == replication_table.5e8$traitpairs[i] &
                                                   !(mr_table.pops$exposure_pop %in% c("UK Biobank", "Biobank Japan")) &
                                                   mr_table.pops$exposure_pop == mr_table.pops$outcome_pop,] %>%
                                                   select("exposure_pop", "grapple5e8.beta", "grapple5e8.se")
  traitpair.estimates.5e8$grapple5e8.p <- 2 * pnorm(abs(traitpair.estimates.5e8$grapple5e8.beta/traitpair.estimates.5e8$grapple5e8.se),
                                                    lower.tail = FALSE)

  if (dim(traitpair.estimates.5e8)[1] >=1 ) {
    rep.pairs.available[i] <- dim(traitpair.estimates.5e8)[1]

    #we grab all summary statistics
    rep1.beta[i] <- traitpair.estimates.5e8$grapple5e8.beta[1]
    rep1.pvalue[i] <- traitpair.estimates.5e8$grapple5e8.p[1]
    rep1.targetpop[i] <- traitpair.estimates.5e8$exposure_pop[1]

    rep2.beta[i] <- NA
    rep2.pvalue[i] <- NA
    rep2.targetpop[i] <- NA

    #we have only 1 or 2 pairs available for replication, which means I can get away with some pretty crappy code
    if (rep.pairs.available[i] == 2) {
        rep2.beta[i] <- traitpair.estimates.5e8$grapple5e8.beta[2]
        rep2.pvalue[i] <- traitpair.estimates.5e8$grapple5e8.p[2]
        rep2.targetpop[i] <- traitpair.estimates.5e8$exposure_pop[2]
    }

    sign_align.5e8 <- traitpair.estimates.5e8[sign(traitpair.estimates.5e8$grapple5e8.beta) == sign(replication_table.5e8$crosspop.5e8[i]),]
    if (dim(sign_align.5e8)[1] >= 1) {
      rep.minp[i] <- min(sign_align.5e8$grapple5e8.p)
    } else { #there are no estimates with the same sign as the cross-population estimate.
      rep.minp[i] <- 1
    }
  } else {
    rep.pairs.available[i] <- 0
    rep1.beta[i] <- rep2.beta[i] <- rep1.pvalue[i] <- rep2.pvalue[i] <- rep1.targetpop[i] <- rep2.targetpop[i] <- NA
    rep.minp[i] <- NA
  }
}

replication_table.5e8$rep.pairs.available <- rep.pairs.available
replication_table.5e8$rep1.beta <- rep1.beta
replication_table.5e8$rep1.pvalue <- rep1.pvalue
replication_table.5e8$rep1.targetpop <- rep1.targetpop
replication_table.5e8$rep2.beta <- rep2.beta
replication_table.5e8$rep2.pvalue <- rep2.pvalue
replication_table.5e8$rep.minp <- rep.minp
replication_table.5e8$rep2.targetpop <- rep2.targetpop

rep.threshold.5e8 <- 0.05/sum(replication_table.5e8$rep.pairs.available > 0)
replication_table.5e8$replicated.5e8 <- ifelse(replication_table.5e8$rep.pairs.available > 0,
                                           ifelse(replication_table.5e8$rep.minp < rep.threshold.5e8, "Replicated", "Not Replicated"),
                                           "Not Tested")


rep_counts.1e5 <- table(replication_table.1e5$replicated.1e5)
rep_table.1e5 <- table(replication_table.1e5$sig_status.1e5, replication_table.1e5$replicated.1e5)

rep_counts.5e8 <- table(replication_table.5e8$replicated.5e8)
rep_table.5e8 <- table(replication_table.5e8$sig_status.5e8, replication_table.5e8$replicated.5e8)

#=Table 2=
table_categories <- c("Instrument Selection Threshold",
                      "Total trait pairs assessed",
                      "Neither estimate significant",
                      "Both estimates significant",
                      "Only Meta-analysis to BBJ estimate significant",
                      "Replication Analyses Performed (Successful Replications), Only Meta-analysis to BBJ significant",
                      "Only BBJ to BBJ estimate significant",
                      "Replication Analyses Performed (Successful Replications), Only BBJ to BBJ significant")

numbers_1e5 <- c("1e-5",
                 sum(!is.na(sakaue_bbj_table$samepop.1e5) & !is.na(sakaue_bbj_table$crosspop.1e5)),
                 sum(sakaue_bbj_table$samepop.p.1e5 >= 0.05/195 & sakaue_bbj_table$crosspop.p.1e5 >= 0.05/195),
                 sum(sakaue_bbj_table$samepop.p.1e5 < 0.05/195 & sakaue_bbj_table$crosspop.p.1e5 < 0.05/195),
                 sum(sakaue_bbj_table$samepop.p.1e5 >= 0.05/195 & sakaue_bbj_table$crosspop.p.1e5 < 0.05/195),
                 paste0(sum(c(rep_table.1e5[1,1], rep_table.1e5[1,3])), " (", rep_table.1e5[1,3], ")"),
                 sum(sakaue_bbj_table$samepop.p.1e5 < 0.05/195 & sakaue_bbj_table$crosspop.p.1e5 >= 0.05/195),
                 paste0(sum(c(rep_table.1e5[2,1], rep_table.1e5[2,3])), " (", rep_table.1e5[2,3], ")"))

numbers_5e8 <- c("5e-8",
                 sum(!is.na(sakaue_bbj_table$samepop.5e8) & !is.na(sakaue_bbj_table$crosspop.5e8)),
                 sum(sakaue_bbj_table$samepop.p.5e8 >= 0.05/195 & sakaue_bbj_table$crosspop.p.5e8 >= 0.05/195),
                 sum(sakaue_bbj_table$samepop.p.5e8 < 0.05/195 & sakaue_bbj_table$crosspop.p.5e8 < 0.05/195),
                 sum(sakaue_bbj_table$samepop.p.5e8 >= 0.05/195 & sakaue_bbj_table$crosspop.p.5e8 < 0.05/195),
                 paste0(sum(c(rep_table.5e8[1,1], rep_table.5e8[1,3])), " (", rep_table.5e8[1,3], ")"),
                 sum(sakaue_bbj_table$samepop.p.5e8 < 0.05/195 & sakaue_bbj_table$crosspop.p.5e8 >= 0.05/195),
                 paste0(sum(c(rep_table.5e8[2,1], rep_table.5e8[2,3])), " (", rep_table.5e8[2,3], ")"))

fwrite(data.frame(table_categories, numbers_1e5, numbers_5e8), col.names = FALSE, "reproducible_figures/table2_sakaue_bbj.csv")

#=Supplemental Table S7=
table_s7 <- replication_table.1e5 %>% select("traitpairs", "sig_status.1e5", "rep.pairs.available",
                                             "rep1.targetpop", "rep2.targetpop", "replicated.1e5")
table_s7$replication.pops <- ifelse(table_s7$rep.pairs.available == 0, "",
                              ifelse(table_s7$rep.pairs.available == 1, ifelse(is.na(table_s7$rep1.targetpop), table_s7$rep2.targetpop,
                                                                                                              table_s7$rep1.targetpop),
                                     paste0(table_s7$rep1.targetpop, ", ", table_s7$rep2.targetpop)))
table_s7$sig_status.1e5 <- ifelse(table_s7$sig_status.1e5 == "Gain of Significance", "Meta analysis Only", "BBJ Only")
table_s7$sig_status.1e5 <- factor(table_s7$sig_status.1e5, levels = c("Meta analysis Only", "BBJ Only")) #sorting for more intuitive ordering
table_s7 <- table_s7[order(table_s7$sig_status.1e5, table_s7$traitpairs, decreasing = FALSE),]
fwrite(table_s7 %>% select("traitpairs", "sig_status.1e5", "replication.pops", "replicated.1e5"),
       "reproducible_figures/tables7_replication1e5.csv")
#=Supplemental Table S8=
table_s8 <- replication_table.5e8 %>% select("traitpairs", "sig_status.5e8", "rep.pairs.available",
                                             "rep1.targetpop", "rep2.targetpop", "replicated.5e8")
table_s8$replication.pops <- ifelse(table_s8$rep.pairs.available == 0, "",
                              ifelse(table_s8$rep.pairs.available == 1, ifelse(is.na(table_s8$rep1.targetpop), table_s8$rep2.targetpop,
                                                                                                              table_s8$rep1.targetpop),
                                     paste0(table_s8$rep1.targetpop, ", ", table_s8$rep2.targetpop)))
table_s8$sig_status.5e8 <- ifelse(table_s8$sig_status.5e8 == "Gain of Significance", "Meta analysis Only", "BBJ Only")
table_s8$sig_status.5e8 <- factor(table_s8$sig_status.5e8, levels = c("Meta analysis Only", "BBJ Only")) #sorting for more intuitive ordering
table_s8 <- table_s8[order(table_s8$sig_status.5e8, table_s8$traitpairs, decreasing = FALSE),]
fwrite(table_s8 %>% select("traitpairs", "sig_status.5e8", "replication.pops", "replicated.5e8"),
       "reproducible_figures/tables8_replication5e8.csv")

#==Results without Steiger filtering for GRAPPLE at both thresholds==

#Assembling table
unfiltered.grapple1e5.beta <- vector()
unfiltered.grapple1e5.se <- vector()

unfiltered.grapple5e8.beta <- vector()
unfiltered.grapple5e8.se <- vector()

for (i in 1:length(new_mr_studypairs)) {
  current.studypair <- new_mr_studypairs[i]

  #==Extracting MR objects==
  unfiltered.grapple1e5.object <- new_mr[[which(names(new_mr) == current.studypair)]]$mr_grapple_object1e5.unfiltered
  unfiltered.grapple5e8.object <- new_mr[[which(names(new_mr) == current.studypair)]]$mr_grapple_object5e8.unfiltered

  #==Extracting beta hat==
  unfiltered.grapple1e5.beta[i] <- ifelse(identical(unfiltered.grapple1e5.object, NA), NA, unfiltered.grapple1e5.object$beta.hat)
  unfiltered.grapple5e8.beta[i] <- ifelse(identical(unfiltered.grapple5e8.object, NA), NA, unfiltered.grapple5e8.object$beta.hat)

  #==Extracting standard errors==
  unfiltered.grapple1e5.se[i] <- ifelse(identical(unfiltered.grapple1e5.object, NA), NA, unfiltered.grapple1e5.object$beta.var %>% sqrt())
  unfiltered.grapple5e8.se[i] <- ifelse(identical(unfiltered.grapple5e8.object, NA), NA, unfiltered.grapple5e8.object$beta.var %>% sqrt())

  if (i %% 1000 == 0) print(i)
}

#=Supplemental Figure 1=
png(file = "reproducible_figures/figures1_unfiltered_vs_steiger.png", width = 1920, height = 1080)
par(mar = c(5.1, 5.1, 4.1, 2.1), mfcol = c(1,2))

plot(unfiltered.grapple1e5.beta, mr_table$grapple1e5.beta,
     main = "MR Estimates before and after Steiger filtering, GRAPPLE (1e-5)", xlab = "MR Estimates without Filtering", ylab = "MR Estimates with Steiger Filtering",
     pch = 16,
     cex.main = 2.4, cex.lab = 2.2, cex.axis = 2)
abline(a = 0, b = 1, lty = 2)
abline(v = 0)
abline(h = 0)

plot(unfiltered.grapple5e8.beta, mr_table$grapple5e8.beta,
     main = "MR Estimates before and after Steiger filtering, GRAPPLE (5e-8)", xlab = "MR Estimates without Filtering", ylab = "MR Estimates with Steiger Filtering",
     pch = 16,
     cex.main = 2.4, cex.lab = 2.2, cex.axis = 2)
abline(a = 0, b = 1, lty = 2)
abline(v = 0)
abline(h = 0)

dev.off()

#=Supplemental Table S6=
mr_unfiltered.table <- data.frame(study_pair = new_mr_studypairs,
                       exposure_id = gsub("(.*)__(.*)", "\\1", new_mr_studypairs),
                       outcome_id = gsub("(.*)__(.*)", "\\2", new_mr_studypairs),
                       grapple1e5.beta = unfiltered.grapple1e5.beta, grapple1e5.se = unfiltered.grapple1e5.se,
                       grapple5e8.beta = unfiltered.grapple5e8.beta, grapple5e8.se = unfiltered.grapple5e8.se,
                       exposure.rescale = mr_table$exposure.rescale, outcome.rescale = mr_table$outcome.rescale)

mr_unfiltered.table$grapple1e5.beta.rescale <- mr_unfiltered.table$grapple1e5.beta/mr_unfiltered.table$exposure.rescale/mr_unfiltered.table$outcome.rescale
mr_unfiltered.table$grapple1e5.se.rescale <- mr_unfiltered.table$grapple1e5.se/mr_unfiltered.table$exposure.rescale/mr_unfiltered.table$outcome.rescale

mr_unfiltered.table$grapple5e8.beta.rescale <- mr_unfiltered.table$grapple5e8.beta/mr_unfiltered.table$exposure.rescale/mr_unfiltered.table$outcome.rescale
mr_unfiltered.table$grapple5e8.se.rescale <- mr_unfiltered.table$grapple5e8.se/mr_unfiltered.table$exposure.rescale/mr_unfiltered.table$outcome.rescale

table_s6 <- mr_table.exactETOS %>% select(study_pair, exposure_id, outcome_id, is.bestmatch, match.type, ETOS)
table_s6 <- inner_join(table_s6, mr_unfiltered.table %>% select(study_pair, grapple1e5.beta.rescale, grapple1e5.se.rescale, grapple5e8.beta.rescale, grapple5e8.se.rescale), by = "study_pair")

#Calculating z-scores of concordance
table_s6$zconcord.grapple1e5 <- NA
table_s6$zconcord.grapple5e8 <- NA

for (i in 1:dim(table_s6)[1]) {
  if (table_s6$is.bestmatch[i]) { #No need to compare the reference estimate to itself!
    next
  }

  #Identifying the reference estimate pair
  ETOS <- table_s6$ETOS[i]
  ETOS.outcome <- table_s6$outcome_id[i]
  ETOS.reference <- table_s6$exposure_id[table_s6$is.bestmatch & table_s6$ETOS == ETOS & table_s6$outcome_id == ETOS.outcome]

  #Reference estimates and SEs for all five MR estimate types, if they exist
  grapple1e5.reference <- table_s6$grapple1e5.beta.rescale[table_s6$outcome_id == ETOS.outcome & table_s6$exposure_id == ETOS.reference]
  grapple1e5.reference.se <- table_s6$grapple1e5.se.rescale[table_s6$outcome_id == ETOS.outcome & table_s6$exposure_id == ETOS.reference]
  grapple5e8.reference <- table_s6$grapple5e8.beta.rescale[table_s6$outcome_id == ETOS.outcome & table_s6$exposure_id == ETOS.reference]
  grapple5e8.reference.se <- table_s6$grapple5e8.se.rescale[table_s6$outcome_id == ETOS.outcome & table_s6$exposure_id == ETOS.reference]

  #The mismatching estimates and SEs for both unfiltered/Steiger
  grapple1e5.mismatch <- table_s6$grapple1e5.beta.rescale[i]
  grapple1e5.mismatch.se <- table_s6$grapple1e5.se.rescale[i]
  grapple5e8.mismatch <- table_s6$grapple5e8.beta.rescale[i]
  grapple5e8.mismatch.se <- table_s6$grapple5e8.se.rescale[i]

  #Calculating concordance Z-scores
  table_s6$zconcord.grapple1e5[i] <- (grapple1e5.mismatch - grapple1e5.reference)/sqrt(grapple1e5.reference.se^2 + grapple1e5.mismatch.se^2)
  table_s6$zconcord.grapple5e8[i] <- (grapple5e8.mismatch - grapple5e8.reference)/sqrt(grapple5e8.reference.se^2 + grapple5e8.mismatch.se^2)
}

#Calculating concordance z-score p-values
table_s6$pconcord.grapple1e5 <- 2 * pnorm(abs(table_s6$zconcord.grapple1e5), lower.tail = FALSE)
table_s6$pconcord.grapple1e5.bh <- p.adjust(table_s6$pconcord.grapple1e5, method = "BH")

table_s6$pconcord.grapple5e8 <- 2 * pnorm(abs(table_s6$zconcord.grapple5e8), lower.tail = FALSE)
table_s6$pconcord.grapple5e8.bh <- p.adjust(table_s6$pconcord.grapple5e8, method = "BH")

#Calculating "mismatch-to-reference" ratio for GRAPPLE estimates

mtr.ratio.1e5 <- vector()
reference.estimate.1e5 <- vector()
reference.se.1e5 <- vector()

mtr.ratio.5e8 <- vector()
reference.estimate.5e8 <- vector()
reference.se.5e8 <- vector()

for (i in 1:dim(table_s6)[1]) {
  if (table_s6$is.bestmatch[i]) { #mismatch/reference ratio only applies to mismatches.
    mtr.ratio.1e5[i] <- NA
    reference.estimate.1e5[i] <- NA
    reference.se.1e5[i] <- NA
    mtr.ratio.5e8[i] <- NA
    reference.estimate.5e8[i] <- NA
    reference.se.5e8[i] <- NA
    next
  }

  mismatch.estimate.1e5 <- table_s6$grapple1e5.beta.rescale[i]
  mismatch.estimate.5e8 <- table_s6$grapple5e8.beta.rescale[i]
  ETOS <- table_s6$ETOS[i]

  reference.estimate.1e5[i] <- table_s6$grapple1e5.beta.rescale[table_s6$is.bestmatch & table_s6$ETOS == ETOS]
  reference.se.1e5[i] <- table_s6$grapple1e5.se.rescale[table_s6$is.bestmatch & table_s6$ETOS == ETOS]
  mtr.ratio.1e5[i] <- mismatch.estimate.1e5/reference.estimate.1e5[i]

  if (!is.na(mismatch.estimate.5e8)) {
    reference.estimate.5e8[i] <- table_s6$grapple5e8.beta.rescale[table_s6$is.bestmatch & table_s6$ETOS == ETOS]
    reference.se.5e8[i] <- table_s6$grapple5e8.se.rescale[table_s6$is.bestmatch & table_s6$ETOS == ETOS]
    mtr.ratio.5e8[i] <- mismatch.estimate.5e8/reference.estimate.5e8[i]
  } else {
    reference.estimate.5e8[i] <- NA
    reference.se.5e8[i] <- NA
    mtr.ratio.5e8[i] <- NA
  }
}

table_categories <- c("Instrument Selection Threshold",
                      "Total Number of Mismatching Estimates with Reference Estimates",
                      "Mismatching Estimates with Non-Significant Concordance z-scores",
                      "Mismatching Estimates with Significant Concordance z-scores, Same Sign as Reference, Smaller Magnitude",
                      "Mismatching Estimates with Significant Concordance z-scores, Same Sign as Reference, Larger Magnitude",
                      "Mismatching Estimates with Significant Concordance z-scores, Opposite Sign as Reference",
                      "Mismatching Estimates with Significant Concordance z-scores, Opposite Sign as Reference, Reference and Mismatch MR p-values < 0.05")

crit.z <- abs(qnorm(0.05/2))
numbers_1e5 <- c("1e-5",
                 sum(!is.na(table_s6$zconcord.grapple1e5)),
                 sum(table_s6$pconcord.grapple1e5.bh >= 0.05, na.rm = TRUE),
                 sum(table_s6$pconcord.grapple1e5.bh < 0.05 & mtr.ratio.1e5 > 0 & mtr.ratio.1e5 < 1, na.rm = TRUE),
                 sum(table_s6$pconcord.grapple1e5.bh < 0.05 & mtr.ratio.1e5 > 1, na.rm = TRUE),
                 sum(table_s6$pconcord.grapple1e5.bh < 0.05 & mtr.ratio.1e5 < 0, na.rm = TRUE),
                 sum(table_s6$pconcord.grapple1e5.bh < 0.05 & mtr.ratio.1e5 < 0 &
                     abs(table_s6$grapple1e5.beta.rescale/table_s6$grapple1e5.se.rescale) > crit.z &
                     abs(reference.estimate.1e5/reference.se.1e5) > crit.z, na.rm = TRUE))

numbers_5e8 <- c("5e-8",
                 sum(!is.na(table_s6$zconcord.grapple5e8)),
                 sum(table_s6$pconcord.grapple5e8.bh >= 0.05, na.rm = TRUE),
                 sum(table_s6$pconcord.grapple5e8.bh < 0.05 & mtr.ratio.5e8 > 0 & mtr.ratio.5e8 < 1, na.rm = TRUE),
                 sum(table_s6$pconcord.grapple5e8.bh < 0.05 & mtr.ratio.5e8 > 1, na.rm = TRUE),
                 sum(table_s6$pconcord.grapple5e8.bh < 0.05 & mtr.ratio.5e8 < 0, na.rm = TRUE),
                 sum(table_s6$pconcord.grapple5e8.bh < 0.05 & mtr.ratio.5e8 < 0 &
                     abs(table_s6$grapple5e8.beta.rescale/table_s6$grapple5e8.se.rescale) > crit.z &
                     abs(reference.estimate.5e8/reference.se.5e8) > crit.z, na.rm = TRUE))

fwrite(data.frame(table_categories, numbers_1e5, numbers_5e8), col.names = FALSE, "reproducible_figures/tables6_zconcord_grapple_unfiltered.csv")
