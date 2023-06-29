set.seed(1)
library(parallel)
library(robustDAA)
outlier <- "outlier1"
ratio <- "ratio1"
# source(paste0("stool_data/", outlier, "/utility.R"))
## define parameters
n_sam <- c(50, 100)
n_taxa <- 500
set_df <- data.frame(n_sam = rep(n_sam, length(n_taxa)), n_taxa = rep(n_taxa, each = length(n_sam)))
signa_den <- c(0.05, 0.2)
set_df <- cbind(apply(set_df, 2, rep, length(signa_den)), rep(signa_den, each = nrow(set_df)))
signa_streng <- 2
set_df <- cbind(apply(set_df, 2, rep, length(signa_streng)), rep(signa_streng, each = nrow(set_df)))
colnames(set_df) <- c("n_sam", "n_taxa", "signa_den", "signa_streng")
set_df <- as.data.frame(set_df)
nset <- nrow(set_df)

################## without confounder ##################
formula <- "u"
for (iter in seq_len(nset)) {
  #### load datasets
  print(iter)
  n_sam <- as.numeric(set_df$n_sam[iter])
  n_taxa <- as.numeric(set_df$n_taxa[iter])
  signa_den <- as.numeric(set_df$signa_den[iter])
  signa_streng <- as.numeric(set_df$signa_streng[iter])
  dta_list <- readRDS(paste0(
    "stool_data/", outlier, "/", ratio, "/datasets/noconf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))

  #### LinDA method
  linda_res <- mclapply(dta_list, function(dta) {
    Y <- dta$otu.tab.sim
    colnames(Y) <- paste0("sample", seq_len(n_sam))
    rownames(Y) <- paste0("taxon", seq_len(n_taxa))
    Z <- dta$covariate
    colnames(Z) <- "u"
    rownames(Z) <- paste0("sample", seq_len(n_sam))
    Z <- as.data.frame(Z)
    Z$u <- as.factor(Z$u)
    res <- LinDA::linda(Y, Z, paste("~", formula))
    rej <- which(res$output[[1]]$reject == TRUE)
    return(rej)
  }, mc.cores = 50)
  saveRDS(linda_res, paste0(
    "stool_data/", outlier, "/", ratio, "/results/linda_noconf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))

  #### LinDA97 method
  linda97_res <- mclapply(dta_list, function(dta) {
    Y <- dta$otu.tab.sim
    colnames(Y) <- paste0("sample", seq_len(n_sam))
    rownames(Y) <- paste0("taxon", seq_len(n_taxa))
    Z <- dta$covariate
    colnames(Z) <- "u"
    rownames(Z) <- paste0("sample", seq_len(n_sam))
    Z <- as.data.frame(Z)
    Z$u <- as.factor(Z$u)
    ## remove samples for Y
    Y_tmp <- winsor.fun(Y, 0.97)
    N_tmp <- colSums(Y_tmp)
    keep_sam <- which(N_tmp >= 1)
    Y <- Y[, keep_sam]
    ## remove samples for Z
    allvars <- colnames(Z)
    Z <- as.data.frame(Z[keep_sam, ])
    colnames(Z) <- allvars
    res <- LinDA::linda(
      otu.tab = Y, meta = Z, formula = paste("~", formula),
      winsor.quan = 0.97
    )
    rej <- which(res$output[[1]]$reject == TRUE)
    return(rej)
  }, mc.cores = 50)
  saveRDS(linda97_res, paste0(
    "stool_data/", outlier, "/", ratio, "/results/linda97_noconf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))

  #### LinDA90 method
  linda90_res <- mclapply(dta_list, function(dta) {
    Y <- dta$otu.tab.sim
    colnames(Y) <- paste0("sample", seq_len(n_sam))
    rownames(Y) <- paste0("taxon", seq_len(n_taxa))
    Z <- dta$covariate
    colnames(Z) <- "u"
    rownames(Z) <- paste0("sample", seq_len(n_sam))
    Z <- as.data.frame(Z)
    Z$u <- as.factor(Z$u)
    ## remove samples for Y
    Y_tmp <- winsor.fun(Y, 0.90)
    N_tmp <- colSums(Y_tmp)
    keep_sam <- which(N_tmp >= 1)
    Y <- Y[, keep_sam]
    ## remove samples for Z
    allvars <- colnames(Z)
    Z <- as.data.frame(Z[keep_sam, ])
    colnames(Z) <- allvars
    res <- LinDA::linda(
      otu.tab = Y, meta = Z, formula = paste("~", formula),
      winsor.quan = 0.90
    )
    rej <- which(res$output[[1]]$reject == TRUE)
    return(rej)
  }, mc.cores = 50)
  saveRDS(linda90_res, paste0(
    "stool_data/", outlier, "/", ratio, "/results/linda90_noconf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))

  #### huber regression
  huber_res <- mclapply(dta_list, function(dta) {
    Y <- dta$otu.tab.sim
    colnames(Y) <- paste0("sample", seq_len(n_sam))
    rownames(Y) <- paste0("taxon", seq_len(n_taxa))
    Z <- dta$covariate
    colnames(Z) <- "u"
    rownames(Z) <- paste0("sample", seq_len(n_sam))
    Z <- as.data.frame(Z)
    Z$u <- as.factor(Z$u)
    res <- rlm_fun(
      otu_tab = Y, meta = Z, formula = paste("~", formula),
      reg_method = "psi.huber"
    )
    rej <- res$index_select
    return(rej)
  }, mc.cores = 50)
  saveRDS(huber_res, paste0(
    "stool_data/", outlier, "/", ratio, "/results/huber_noconf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))

  #### bisquare regression
  bisquare_res <- mclapply(dta_list, function(dta) {
    Y <- dta$otu.tab.sim
    colnames(Y) <- paste0("sample", seq_len(n_sam))
    rownames(Y) <- paste0("taxon", seq_len(n_taxa))
    Z <- dta$covariate
    colnames(Z) <- "u"
    rownames(Z) <- paste0("sample", seq_len(n_sam))
    Z <- as.data.frame(Z)
    Z$u <- as.factor(Z$u)
    res <- rlm_fun(
      otu_tab = Y, meta = Z, formula = paste("~", formula),
      reg_method = "psi.bisquare"
    )
    rej <- res$index_select
    return(rej)
  }, mc.cores = 50)
  saveRDS(bisquare_res, paste0(
    "stool_data/", outlier, "/", ratio, "/results/bisquare_noconf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))

  #### quantile regression
  qr_res <- mclapply(dta_list, function(dta) {
    Y <- dta$otu.tab.sim
    colnames(Y) <- paste0("sample", seq_len(n_sam))
    rownames(Y) <- paste0("taxon", seq_len(n_taxa))
    Z <- dta$covariate
    colnames(Z) <- "u"
    rownames(Z) <- paste0("sample", seq_len(n_sam))
    Z <- as.data.frame(Z)
    Z$u <- as.factor(Z$u)
    res <- qr_fun(otu_tab = Y, meta = Z, paste("~", formula))
    rej <- res$index_select
    return(rej)
  }, mc.cores = 50)
  saveRDS(qr_res, paste0(
    "stool_data/", outlier, "/", ratio, "/results/qr_noconf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))
}

################## with confounder ##################
formula <- "u + z1 + z2"
for (iter in seq_len(nset)) {
  #### load datasets
  print(iter)
  n_sam <- as.numeric(set_df$n_sam[iter])
  n_taxa <- as.numeric(set_df$n_taxa[iter])
  signa_den <- as.numeric(set_df$signa_den[iter])
  signa_streng <- as.numeric(set_df$signa_streng[iter])
  dta_list <- readRDS(paste0(
    "stool_data/", outlier, "/", ratio, "/datasets/conf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))

  #### LinDA method
  linda_res <- mclapply(dta_list, function(dta) {
    Y <- dta$otu.tab.sim
    colnames(Y) <- paste0("sample", seq_len(n_sam))
    rownames(Y) <- paste0("taxon", seq_len(n_taxa))
    Z <- cbind(dta$covariate, dta$confounder)
    colnames(Z) <- c("u", "z1", "z2")
    rownames(Z) <- paste0("sample", seq_len(n_sam))
    Z <- as.data.frame(Z)
    Z$u <- as.factor(Z$u)
    Z$z2 <- as.factor(Z$z2)
    res <- LinDA::linda(Y, Z, paste("~", formula))
    rej <- which(res$output[[1]]$reject == TRUE)
    return(rej)
  }, mc.cores = 50)
  saveRDS(linda_res, paste0(
    "stool_data/", outlier, "/", ratio, "/results/linda_conf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))

  #### LinDA97 method
  linda97_res <- mclapply(dta_list, function(dta) {
    Y <- dta$otu.tab.sim
    colnames(Y) <- paste0("sample", seq_len(n_sam))
    rownames(Y) <- paste0("taxon", seq_len(n_taxa))
    Z <- cbind(dta$covariate, dta$confounder)
    colnames(Z) <- c("u", "z1", "z2")
    rownames(Z) <- paste0("sample", seq_len(n_sam))
    Z <- as.data.frame(Z)
    Z$u <- as.factor(Z$u)
    Z$z2 <- as.factor(Z$z2)
    ## remove samples for Y
    Y_tmp <- winsor.fun(Y, 0.97)
    N_tmp <- colSums(Y_tmp)
    keep_sam <- which(N_tmp >= 1)
    Y <- Y[, keep_sam]
    ## remove samples for Z
    allvars <- colnames(Z)
    Z <- as.data.frame(Z[keep_sam, ])
    colnames(Z) <- allvars
    res <- LinDA::linda(
      otu.tab = Y, meta = Z, formula = paste("~", formula),
      winsor.quan = 0.97
    )
    rej <- which(res$output[[1]]$reject == TRUE)
    return(rej)
  }, mc.cores = 50)
  saveRDS(linda97_res, paste0(
    "stool_data/", outlier, "/", ratio, "/results/linda97_conf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))

  #### LinDA90 method
  linda90_res <- mclapply(dta_list, function(dta) {
    Y <- dta$otu.tab.sim
    colnames(Y) <- paste0("sample", seq_len(n_sam))
    rownames(Y) <- paste0("taxon", seq_len(n_taxa))
    Z <- cbind(dta$covariate, dta$confounder)
    colnames(Z) <- c("u", "z1", "z2")
    rownames(Z) <- paste0("sample", seq_len(n_sam))
    Z <- as.data.frame(Z)
    Z$u <- as.factor(Z$u)
    Z$z2 <- as.factor(Z$z2)
    ## remove samples for Y
    Y_tmp <- winsor.fun(Y, 0.90)
    N_tmp <- colSums(Y_tmp)
    keep_sam <- which(N_tmp >= 1)
    Y <- Y[, keep_sam]
    ## remove samples for Z
    allvars <- colnames(Z)
    Z <- as.data.frame(Z[keep_sam, ])
    colnames(Z) <- allvars
    res <- LinDA::linda(
      otu.tab = Y, meta = Z, formula = paste("~", formula),
      winsor.quan = 0.90
    )
    rej <- which(res$output[[1]]$reject == TRUE)
    return(rej)
  }, mc.cores = 50)
  saveRDS(linda90_res, paste0(
    "stool_data/", outlier, "/", ratio, "/results/linda90_conf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))

  #### huber regression
  huber_res <- mclapply(dta_list, function(dta) {
    Y <- dta$otu.tab.sim
    colnames(Y) <- paste0("sample", seq_len(n_sam))
    rownames(Y) <- paste0("taxon", seq_len(n_taxa))
    Z <- cbind(dta$covariate, dta$confounder)
    colnames(Z) <- c("u", "z1", "z2")
    rownames(Z) <- paste0("sample", seq_len(n_sam))
    Z <- as.data.frame(Z)
    Z$u <- as.factor(Z$u)
    Z$z2 <- as.factor(Z$z2)
    res <- rlm_fun(
      otu_tab = Y, meta = Z, formula = paste("~", formula),
      reg_method = "psi.huber"
    )
    rej <- res$index_select
    return(rej)
  }, mc.cores = 50)
  saveRDS(huber_res, paste0(
    "stool_data/", outlier, "/", ratio, "/results/huber_conf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))

  #### bisquare regression
  bisquare_res <- mclapply(dta_list, function(dta) {
    Y <- dta$otu.tab.sim
    colnames(Y) <- paste0("sample", seq_len(n_sam))
    rownames(Y) <- paste0("taxon", seq_len(n_taxa))
    Z <- cbind(dta$covariate, dta$confounder)
    colnames(Z) <- c("u", "z1", "z2")
    rownames(Z) <- paste0("sample", seq_len(n_sam))
    Z <- as.data.frame(Z)
    Z$u <- as.factor(Z$u)
    Z$z2 <- as.factor(Z$z2)
    res <- rlm_fun(
      otu_tab = Y, meta = Z, formula = paste("~", formula),
      reg_method = "psi.bisquare"
    )
    rej <- res$index_select
    return(rej)
  }, mc.cores = 50)
  saveRDS(bisquare_res, paste0(
    "stool_data/", outlier, "/", ratio, "/results/bisquare_conf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))

  #### quantile regression
  qr_res <- mclapply(dta_list, function(dta) {
    Y <- dta$otu.tab.sim
    colnames(Y) <- paste0("sample", seq_len(n_sam))
    rownames(Y) <- paste0("taxon", seq_len(n_taxa))
    Z <- cbind(dta$covariate, dta$confounder)
    colnames(Z) <- c("u", "z1", "z2")
    rownames(Z) <- paste0("sample", seq_len(n_sam))
    Z <- as.data.frame(Z)
    Z$u <- as.factor(Z$u)
    Z$z2 <- as.factor(Z$z2)
    res <- qr_fun(otu_tab = Y, meta = Z, paste("~", formula))
    rej <- res$index_select
    return(rej)
  }, mc.cores = 50)
  saveRDS(qr_res, paste0(
    "stool_data/", outlier, "/", ratio, "/results/qr_conf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))
}
