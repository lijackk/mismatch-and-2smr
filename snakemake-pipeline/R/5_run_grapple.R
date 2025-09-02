library(dplyr)
library(purrr)
library(GRAPPLE)
library(stringr)


out <- snakemake@output[["out"]]
R_est_file <- snakemake@input[["R"]]
pt <- as.numeric(snakemake@wildcards[["pt"]])
z_files = unlist(snakemake@input[["Z"]])
#niter = snakemake@params[["max_iter"]]
#seed <- snakemake@wildcards[["fs"]]
#set.seed(seed)
niter = 20


# Read in data
X <- map_dfr(z_files, readRDS)

ntrait <- X %>%
  dplyr::select(ends_with(".z")) %>%
  ncol()

Z_hat <- X %>%
  dplyr::select(ends_with(".z")) %>%
  as.matrix()

se_hat <- X %>%
  dplyr::select(ends_with(".se")) %>%
  as.matrix()

beta_hat <- Z_hat*se_hat
pvals <- 2*pnorm(-abs(Z_hat))

nms <- names(X)[grep(".z$", names(X))] %>% str_replace(".z$", "")
nms2 <- names(X)[grep(".se$", names(X))] %>% str_replace(".se$", "")
stopifnot(all(nms == nms2))

R <- readRDS(R_est_file)
stopifnot(all(R$names  == nms))
# z_order <- match(R$names, nms)
# se_hat <- se_hat[,z_order]
# beta_hat <- beta_hat[,z_order]

Rcor <- cov2cor(R$R)
p <- ncol(beta_hat)

grapple_dat <- data.frame(cbind(X$snp, beta_hat, se_hat))
names(grapple_dat) <- c("snp", "gamma_out", paste0("gamma_exp", 1:(p-1)),
                        "se_out", paste0("se_exp", 1:(p-1)))

grapple_dat[,2:dim(grapple_dat)[2]] <- apply(grapple_dat[,2:dim(grapple_dat)[2]], 2, function(x){as.numeric(x)})

grapple_dat$selection_pvals <- apply(pvals[,-1, drop = F],1, min)

ix <- which(grapple_dat$selection_pvals < pt)

if (length(ix) == 0) print("No sufficiently strong instruments were found!")

if (length(ix) > 0) { #UPDATE: 7/7/2024: this section is only run if ANY strong-enough instruments were found.
  if(length(ix) > 5000){
     ixn <- sample(ix, size = 5000, replace = FALSE)
     grapple_dat <- grapple_dat[ixn,]
  }else{
     grapple_dat <- grapple_dat[ix,]
  }
  
  #UPDATE, 6/30/2024: if grapple_dat is 7 SNPs or fewer, we can't run residual diagnostics; so I've created an if/else statement for that instead.
  
  if (dim(grapple_dat)[1] > 7) {
    t <- system.time(
    res <- grappleRobustEst(data = grapple_dat,
                            plot.it =FALSE,
                            p.thres = pt,
                            cor.mat = Rcor,
                            loss.function = "tukey",
                            niter = 20)) #UPDATE: 7/4/2024 - 1 million max iterations? Really?
  } else {
    t <- system.time(
    res <- grappleRobustEst(data = grapple_dat,
                            plot.it =FALSE,
                            p.thres = pt,
                            cor.mat = Rcor,
                            niter = 20,
                            loss.function = "tukey",
                            diagnosis = FALSE)) #no residual diagnostics.
  }
  
  res$time <- t
  res$names <- nms
  
  saveRDS(list(grapple_data = grapple_dat, result = res), file=out)
}