set.seed(1)
library(parallel)
outlier <- "outlier0"
source(paste0("loglinear_data/", outlier, "/utility.R"))
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

  #### LinDA method
  linda_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
    res <- LinDA::linda(Y, Z, paste("~", formula))
    rej <- which(res$output[[1]]$reject == TRUE)
    return(rej)
  }, mc.cores = 50)
  ## save results
  saveRDS(linda_res, paste0(
    "loglinear_data/", outlier, "/results/linda_nocon_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))

  #### LinDA97 method
  linda97_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
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
    "loglinear_data/", outlier, "/results/linda97_nocon_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))

  #### LinDA90 method
  linda90_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
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
    "loglinear_data/", outlier, "/results/linda90_nocon_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))

  # #### LinDA winsor method
  # linda_winsor_res <- mclapply(dta_list, function(dta) {
  #   Y <- dta$Y
  #   Z <- dta$Z
  #   res <- linda_winsor(Y, Z, paste("~", formula))
  #   rej <- res$index_select
  #   return(rej)
  # }, mc.cores = 50)
  # ## save results
  # saveRDS(linda_winsor_res, paste0(
  #   "loglinear_data/", outlier, "/results/linda_winsor_nocon_n", n,
  #   "gamma", gamma, "mu", mu_use, ".rds"
  # ))

  #### Huber method
  huber_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
    res <- rlm_fun(
      Y = Y, Z = Z, formula = paste("~", formula),
      res_method = "psi.huber",
      test_method = "t", adj_method = "BH"
    )
    rej <- res$index_select
    return(rej)
  }, mc.cores = 50)
  ## save results
  saveRDS(huber_res, paste0(
    "loglinear_data/", outlier, "/results/huber_nocon_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))

  #### bisquare method
  bisquare_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
    res <- rlm_fun(
      Y = Y, Z = Z, formula = paste("~", formula),
      res_method = "psi.bisquare",
      test_method = "t", adj_method = "BH"
    )
    rej <- res$index_select
    return(rej)
  }, mc.cores = 50)
  ## save results
  saveRDS(bisquare_res, paste0(
    "loglinear_data/", outlier, "/results/bisquare_nocon_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))

  #### quantile regression
  qr_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
    res <- qr_fun(Y, Z, paste("~", formula))
    rej <- res$index_select
    return(rej)
  }, mc.cores = 50)
  ## save results
  saveRDS(qr_res, paste0(
    "loglinear_data/", outlier, "/results/qr_nocon_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))
}

################## with confounder ##################
formula <- "u + z1 + z2"
for (iter_para in seq_len(n_setting)) {
  print(iter_para)
  #### load datasets
  para <- setting[iter_para, ]
  n <- para[1]
  gamma <- para[2]
  mu_use <- para[3]
  dta_list <- readRDS(paste0(
    "loglinear_data/", outlier, "/datasets/con_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))

  #### LinDA method
  linda_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
    res <- LinDA::linda(Y, Z, paste("~", formula))
    rej <- which(res$output[[1]]$reject == TRUE)
    return(rej)
  }, mc.cores = 50)
  ## save results
  saveRDS(linda_res, paste0(
    "loglinear_data/", outlier, "/results/linda_con_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))

  #### LinDA97 method
  linda97_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
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
    "loglinear_data/", outlier, "/results/linda97_con_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))

  #### LinDA90 method
  linda90_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
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
    "loglinear_data/", outlier, "/results/linda90_con_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))

  # #### LinDA winsor method
  # linda_winsor_res <- mclapply(dta_list, function(dta) {
  #   Y <- dta$Y
  #   Z <- dta$Z
  #   res <- linda_winsor(Y, Z, paste("~", formula))
  #   rej <- res$index_select
  #   return(rej)
  # }, mc.cores = 50)
  # ## save results
  # saveRDS(linda_winsor_res, paste0(
  #   "loglinear_data/", outlier, "/results/linda_winsor_con_n", n,
  #   "gamma", gamma, "mu", mu_use, ".rds"
  # ))

  #### Huber method
  huber_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
    res <- rlm_fun(
      Y = Y, Z = Z, formula = paste("~", formula),
      res_method = "psi.huber",
      test_method = "t", adj_method = "BH"
    )
    rej <- res$index_select
    return(rej)
  }, mc.cores = 50)
  ## save results
  saveRDS(huber_res, paste0(
    "loglinear_data/", outlier, "/results/huber_con_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))

  #### bisquare method
  bisquare_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
    res <- rlm_fun(
      Y = Y, Z = Z, formula = paste("~", formula),
      res_method = "psi.bisquare",
      test_method = "t", adj_method = "BH"
    )
    rej <- res$index_select
    return(rej)
  }, mc.cores = 50)
  ## save results
  saveRDS(bisquare_res, paste0(
    "loglinear_data/", outlier, "/results/bisquare_con_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))

  #### quantile regression
  qr_res <- mclapply(dta_list, function(dta) {
    Y <- dta$Y
    Z <- dta$Z
    res <- qr_fun(Y, Z, paste("~", formula))
    rej <- res$index_select
    return(rej)
  }, mc.cores = 50)
  ## save results
  saveRDS(qr_res, paste0(
    "loglinear_data/", outlier, "/results/qr_con_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))
}
