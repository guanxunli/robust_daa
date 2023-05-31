set.seed(1)
library(parallel)
outlier <- "outlier2"
source(paste0("mixeffect_data/", outlier, "/utility.R"))
## load parameters
para0 <- readRDS(paste0("mixeffect_data/", outlier, "/datasets/log.normal.para.rds"))
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
mix_vec <- c("mix1", "mix2")
setting <- cbind(apply(setting, 2, rep, length(mix_vec)), rep(mix_vec, each = nrow(setting)))
colnames(setting) <- c("sample.size", "sig.density", "sig.strength", "model")
n_setting <- nrow(setting)

################## simulations ##################
formula <- "u + (1|id)"
for (iter_para in seq_len(n_setting)) {
  #### load datasets
  print(iter_para)
  para <- setting[iter_para, ]
  n <- as.numeric(para[1])
  gamma <- as.numeric(para[2])
  mu_use <- as.numeric(para[3])
  model <- para[4]
  dta_list <- readRDS(paste0(
    "mixeffect_data/", outlier, "/datasets/n", n,
    "gamma", gamma, "mu", mu_use, model, ".rds"
  ))

  #### LinDA method
  linda_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
    Z$id <- as.factor(Z$id)
    res <- LinDA::linda(Y, Z, paste("~", formula))
    rej <- which(res$output[[1]]$reject == TRUE)
    return(rej)
  }, mc.cores = 50)
  ## save results
  saveRDS(linda_res, paste0(
    "mixeffect_data/", outlier, "/results/linda_n", n,
    "gamma", gamma, "mu", mu_use, model, ".rds"
  ))

  #### LinDA 97
  linda97_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
    Z$id <- as.factor(Z$id)
    ## remove samples for Y
    Y_tmp <- winsor.fun(Y, 0.97)
    N_tmp <- colSums(Y_tmp)
    keep_sam <- which(N_tmp >= 1)
    Y <- Y[, keep_sam]
    ## remove samples for Z
    allvars <- colnames(Z)
    Z <- as.data.frame(Z[keep_sam, ])
    colnames(Z) <- allvars
    ## linda method
    res <- LinDA::linda(
      otu.tab = Y, meta = Z, formula = paste("~", formula),
      winsor.quan = 0.97
    )
    rej <- which(res$output[[1]]$reject == TRUE)
    return(rej)
  }, mc.cores = 50)
  ## save results
  saveRDS(linda97_res, paste0(
    "mixeffect_data/", outlier, "/results/linda97_n", n,
    "gamma", gamma, "mu", mu_use, model, ".rds"
  ))

  #### LinDA 90
  linda90_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
    Z$id <- as.factor(Z$id)
    ## remove samples for Y
    Y_tmp <- winsor.fun(Y, 0.90)
    N_tmp <- colSums(Y_tmp)
    keep_sam <- which(N_tmp >= 1)
    Y <- Y[, keep_sam]
    ## remove samples for Z
    allvars <- colnames(Z)
    Z <- as.data.frame(Z[keep_sam, ])
    colnames(Z) <- allvars
    ## linda method
    res <- LinDA::linda(
      otu.tab = Y, meta = Z, formula = paste("~", formula),
      winsor.quan = 0.90
    )
    rej <- which(res$output[[1]]$reject == TRUE)
    return(rej)
  }, mc.cores = 50)
  ## save results
  saveRDS(linda90_res, paste0(
    "mixeffect_data/", outlier, "/results/linda90_n", n,
    "gamma", gamma, "mu", mu_use, model, ".rds"
  ))

  #### Huber method
  huber_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
    Z$id <- as.factor(Z$id)
    res <- rlm_fun(Y, Z, paste("~", formula))
    rej <- res$index_select
    return(rej)
  }, mc.cores = 50)
  ## save results
  saveRDS(huber_res, paste0(
    "mixeffect_data/", outlier, "/results/huber_n", n,
    "gamma", gamma, "mu", mu_use, model, ".rds"
  ))
}
