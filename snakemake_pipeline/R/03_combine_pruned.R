## Combining the pruned datasets into one full dataset, estimating residual correlation R, then 

source("renv/activate.R")
library(dplyr)
library(purrr)
library(GRAPPLE)
library(stringr)
library(GFA)

sumstat_output_path <- snakemake@output[["sumstat_final"]]
Restimate_output_path <- snakemake@output[["Restimate"]]

sumstat_pthresh <- 1e-5
z_files = unlist(snakemake@input[["sumstat_prune"]])
cond_num = 1e3
mwc <- TRUE

##===Reading in the z-score files===
X <- map_dfr(z_files, readRDS)

##===Extracting individual columns and the names of the dataset===

Z_hat <- X %>%
  select(ends_with(".z")) %>%
  as.matrix()

se_hat <- X %>%
  dplyr::select(ends_with(".se")) %>%
  as.matrix()

beta_hat <- Z_hat*se_hat
pvals <- 2*pnorm(-abs(Z_hat))

af <- X %>%
  dplyr::select(ends_with(".af"))

nms <- colnames(Z_hat) %>% stringr::str_replace(".z$", "")

#===Calculating residual correlation via p-value thresholding===
Rpt <- R_pt(B_hat = Z_hat,
            S_hat = matrix(1, nrow = nrow(Z_hat), ncol = ncol(Z_hat)),
            p_val_thresh = 0.05,
            return_cor = TRUE,
            make_well_conditioned = mwc,
            cond_num = cond_num
            )

ret <- list(R = Rpt, names = nms)
saveRDS(ret, file = Restimate_output_path)

#===Keeping only the variants with a p-value of at least 1e-5===
sumstat_final <- data.frame(cbind(X$snp, beta_hat, se_hat, pvals, af))
names(sumstat_final) <- c("SNP", "bx", "by", "bxse", "byse", "pval_x", "pval_y", "af_x", "af_y")

sumstat_final$selection_pvals <- apply(pvals[,-2, drop = F],1, min) ##the outcome is now the second entry, not the first
ix <- which(sumstat_final$selection_pvals < 1e-5)

if (length(ix) == 0) {
  print("No sufficiently strong instruments were found!")
  sumstat_final <- sumstat_final[ix,]
}

if (length(ix) > 0) {
  if(length(ix) > 5000){
     ixn <- sample(ix, size = 5000, replace = FALSE)
     sumstat_final <- sumstat_final[ixn,]
  }else{
     sumstat_final <- sumstat_final[ix,]
  }
}
saveRDS(sumstat_final, file = sumstat_output_path)  