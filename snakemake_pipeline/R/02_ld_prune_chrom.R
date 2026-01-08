source("renv/activate.R")
library(dplyr)
library(ieugwasr)

X <- readRDS(snakemake@input[["sumstat_combined"]])
print(X)
r2_thresh <- 0.01
clump_kb <- 1000
ref_path <- snakemake@input[["LD"]]
print(ref_path)
bfile <- gsub("\\.bim","", ref_path[1])
print(bfile)
out <- snakemake@output[["sumstat_prune"]]
print(out)
pthresh <- 1

Z_hat <- X %>%
  select(ends_with(".z")) %>%
  as.matrix()
Z_hat <- Z_hat[,-2,drop = FALSE] ##the outcome is now the second entry, not the first

zmax <- apply(abs(Z_hat), 1, function(x){max(x, na.rm=T)})
myp <- 2*pnorm(-abs(zmax))

dat <- data.frame(rsid = X$snp, pval = myp)

dat_clump <- ld_clump(dat = dat,
                     clump_r2 = r2_thresh,
                     clump_p = pthresh,
                     clump_kb = clump_kb,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile)

ix <- which(X$snp %in% dat_clump$rsid)
X <- X[ix,]

saveRDS(X, file=out)