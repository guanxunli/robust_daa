set.seed(1)
library(parallel)
outlier <- "outlier3"
source("other_methods/utility.R")
## load parameters
para0 <- readRDS(paste0("loglinear_data/", outlier, "/datasets/log.normal.para.rds"))
beta0 <- para0$beta0
sigma2 <- para0$sigma2
# parameter use
sample.size.vec <- c(50, 200)
m <- 500
n_sim <- 100
# define settings
sig.density.vec <- c(0.05, 0.2)
sig.strength.vec <- seq(1.05, 2, length.out = 6)
s1 <- 2
s2 <- 2
s3 <- 6
sample.size <- rep(sample.size.vec, each = s2 * s3)
sig.density <- rep(rep(sig.density.vec, each = s3), s1)
sig.strength <- rep(sig.strength.vec, s1 * s2)
setting <- cbind(sample.size, sig.density, sig.strength)
n_setting <- nrow(setting)

################## without confounder ##################
formula <- "u"
for (iter_para in seq_len(n_setting)) {
  print(iter_para)
  #### load datasets
  para <- setting[iter_para, ]
  n <- para[1]
  gamma <- para[2]
  mu_use <- para[3]
  dta_list <- readRDS(paste0(
    "loglinear_data/", outlier, "/datasets/nocon_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))
  
  #### aldex2 method
  aldex2_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
    rej <- try(aldex2.fun(Y, Z, formula, alpha = 0.05), silent = TRUE)
    if (class(rej) == "try-error") {
      rej <- NULL
    }
    return(rej)
  }, mc.cores = 50)
  ## save results
  saveRDS(aldex2_res, paste0(
    "other_methods/", outlier, "/results/aldex2_nocon_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))
  
  #### ancombc method
  ancombc_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
    rej <- try(ancombc.fun(Y, Z, formula, alpha = 0.05), silent = TRUE)
    if (class(rej) == "try-error") {
      rej <- NULL
    }
    return(rej)
  }, mc.cores = 50)
  ## save results
  saveRDS(ancombc_res, paste0(
    "other_methods/", outlier, "/results/ancombc_nocon_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))
  
  #### maaslin2 method
  maaslin2_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
    rej <- try(maaslin2.fun(Y, Z, formula, alpha = 0.05), silent = TRUE)
    if (class(rej) == "try-error") {
      rej <- NULL
    }
    return(rej)
  }, mc.cores = 50)
  ## save results
  saveRDS(maaslin2_res, paste0(
    "other_methods/", outlier, "/results/maaslin2_nocon_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))
}
