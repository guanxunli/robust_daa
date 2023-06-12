####################### robust regression #######################
rlm_fun <- function(Y, Z, formula, adaptive = TRUE, imputation = FALSE,
                    pseudo_cnt = 0.5, corr_cut = 0.1, prev_cut = 0, lib_cut = 1,
                    res_method = "psi.huber", n_cv = 5, hyper_para_vec = NULL,
                    test_method = "t", adj_method = "BH", alpha = 0.05) {
  #### check input
  if (any(is.na(Y))) {
    stop("The OTU table contains NAs! Please remove!\n")
  }
  allvars <- colnames(Z)
  
  ## pre processing
  keep.sam <- which(colSums(Y) >= lib_cut & rowSums(is.na(Z)) == 0)
  Y <- Y[, keep.sam]
  Z <- as.data.frame(Z[keep.sam, ])
  colnames(Z) <- allvars
  
  n <- ncol(Y)
  keep.tax <- which(rowSums(Y > 0) / n >= prev_cut)
  Y <- Y[keep.tax, ]
  m <- nrow(Y)
  
  ## some samples may have zero total counts after screening taxa
  if (any(colSums(Y) == 0)) {
    ind <- which(colSums(Y) > 0)
    Y <- Y[, ind]
    Z <- as.data.frame(Z[ind, ])
    names(Z) <- allvars
    keep.sam <- keep.sam[ind]
    n <- ncol(Y)
  }
  
  ## scaling numerical variables
  ind <- sapply(seq_len(ncol(Z)), function(i) is.numeric(Z[, i]))
  Z[, ind] <- scale(Z[, ind])
  p <- ncol(Z) + 1
  
  ## handling zeros
  if (any(Y == 0)) {
    N <- colSums(Y)
    if (adaptive) {
      logN <- log(N)
      tmp <- lm(as.formula(paste0("logN", formula)), Z)
      corr.pval <- coef(summary(tmp))[-1, "Pr(>|t|)"]
      if (any(corr.pval <= corr_cut)) {
        imputation <- TRUE
      } else {
        imputation <- FALSE
      }
    }
    if (imputation) {
      N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
      N.mat[Y > 0] <- 0
      tmp <- N[max.col(N.mat)]
      Y <- Y + N.mat / tmp
    } else {
      Y <- Y + pseudo_cnt
    }
  }
  
  ## CLR transformation
  logY <- log2(Y)
  W <- t(logY) - colMeans(logY)
  
  ## robust linear regression
  options(warn = -1)
  formula_use <- as.formula(paste0("w_tmp", formula))
  formula_train <- as.formula(paste0("w_train", formula))
  alpha_vec <- numeric(m)
  sd_alpha <- numeric(m)
  
  ## Huber method
  if (res_method == "psi.huber") {
    if (is.null(hyper_para_vec)) {
      hyper_para_vec <- exp(seq(log(0.1), log(1.345), length = 10))
    }
    n_hyper <- length(hyper_para_vec)
    cv_fold <- origami::folds_vfold(n = n, V = n_cv)
    for (iter_m in seq_len(m)) {
      ## sample index
      w_tmp <- W[, iter_m]
      resid_para <- numeric(n_hyper)
      resid_cv <- numeric(n_cv)
      for (iter_hyper in seq_len(n_hyper)) {
        hyper_para_use <- hyper_para_vec[iter_hyper]
        for (iter_fold in seq_len(n_cv)) {
          ## load index
          cv_tmp <- cv_fold[[iter_fold]]
          index_train <- cv_tmp$training_set
          index_cross <- cv_tmp$validation_set
          ## load data
          Z_cross <- data.frame(Z[index_cross, ])
          colnames(Z_cross) <- colnames(Z)
          w_cross <- w_tmp[index_cross]
          w_train <- w_tmp[index_train]
          Z_train <- data.frame(Z[index_train, ])
          colnames(Z_train) <- colnames(Z)
          ## fit model
          fit_tmp <- MASS::rlm(formula_train,
                               data = Z_train, psi = res_method, maxit = 50,
                               k = hyper_para_use
          )
          resid_tmp <- w_cross - predict(fit_tmp, data = Z_cross)
          resid_cv[iter_fold] <- sum(resid_tmp * resid_tmp)
        }
        resid_para[iter_hyper] <- mean(resid_cv)
      }
      ## fit model
      hyper_para <- hyper_para_vec[which.min(resid_para)]
      fit <- MASS::rlm(formula_use, data = Z, psi = res_method, maxit = 50, k = hyper_para)
      alpha_vec[iter_m] <- coef(fit)["u1"]
      summary_fit <- summary(fit)
      sd_alpha[iter_m] <- coef(summary_fit)["u1", 2]
    }
  } else if (res_method == "psi.bisquare") {
    ## bisquare method
    if (is.null(hyper_para_vec)) {
      hyper_para_vec <- exp(seq(log(1), log(4.685), length = 10))
    }
    n_hyper <- length(hyper_para_vec)
    cv_fold <- origami::folds_vfold(n = n, V = n_cv)
    for (iter_m in seq_len(m)) {
      ## sample index
      w_tmp <- W[, iter_m]
      resid_para <- numeric(n_hyper)
      resid_cv <- numeric(n_cv)
      for (iter_hyper in seq_len(n_hyper)) {
        hyper_para_use <- hyper_para_vec[iter_hyper]
        for (iter_fold in seq_len(n_cv)) {
          ## load index
          cv_tmp <- cv_fold[[iter_fold]]
          index_train <- cv_tmp$training_set
          index_cross <- cv_tmp$validation_set
          ## load data
          Z_cross <- data.frame(Z[index_cross, ])
          colnames(Z_cross) <- colnames(Z)
          w_cross <- w_tmp[index_cross]
          w_train <- w_tmp[index_train]
          Z_train <- data.frame(Z[index_train, ])
          colnames(Z_train) <- colnames(Z)
          ## fit model
          fit_tmp <- MASS::rlm(formula_train,
                               data = Z_train, psi = res_method, maxit = 50,
                               c = hyper_para_use
          )
          resid_tmp <- w_cross - predict(fit_tmp, data = Z_cross)
          resid_cv[iter_fold] <- sum(resid_tmp * resid_tmp)
        }
        resid_para[iter_hyper] <- mean(resid_cv)
      }
      ## fit model
      hyper_para <- hyper_para_vec[which.min(resid_para)]
      fit <- MASS::rlm(formula_use, data = Z, psi = res_method, maxit = 50, c = hyper_para)
      alpha_vec[iter_m] <- coef(fit)["u1"]
      summary_fit <- summary(fit)
      sd_alpha[iter_m] <- coef(summary_fit)["u1", 2]
    }
  }
  options(warn = 0)
  
  #### mode correction
  # bias <- median(alpha_vec)
  bias <- modeest::mlv(sqrt(n) * alpha_vec,
                       method = "meanshift", kernel = "gaussian"
  ) / sqrt(n)
  alpha_correct <- alpha_vec - bias
  
  #### test
  if (test_method == "t") {
    ## t-test
    stat <- alpha_correct / sd_alpha
    p_value <- 2 * pt(-abs(stat), df = n - p)
    p_adj <- p.adjust(p_value, method = adj_method)
    index_select <- which(p_adj < alpha)
  } else if (test_method == "chi") {
    ## chi square test
    stat <- alpha_correct * alpha_correct / (sd_alpha * sd_alpha)
    p_value <- 1 - pchisq(stat, df = 1)
    p_adj <- p.adjust(p_value, method = adj_method)
    index_select <- which(p_adj < alpha)
  }
  
  #### return results
  return(list(
    alpha_correct = alpha_correct, sd_alpha = sd_alpha, p_value = p_value,
    index_select = index_select
  ))
}

####################### quantile regression #######################
qr_fun <- function(Y, Z, formula, tau_vec = NULL, adaptive = TRUE, imputation = FALSE,
                   pseudo_cnt = 0.5, corr_cut = 0.1, prev_cut = 0, lib_cut = 1, n_cv = 5,
                   test_method = "t", adj_method = "BH", alpha = 0.05) {
  #### check input
  if (any(is.na(Y))) {
    stop("The OTU table contains NAs! Please remove!\n")
  }
  allvars <- colnames(Z)
  if (is.null(tau_vec)) {
    tau_vec <- seq(0.2, 0.8, by = 0.1)
    n_tau <- length(tau_vec)
  } else {
    n_tau <- length(tau_vec)
  }
  
  ## preprocessing
  keep.sam <- which(colSums(Y) >= lib_cut & rowSums(is.na(Z)) == 0)
  Y <- Y[, keep.sam]
  Z <- as.data.frame(Z[keep.sam, ])
  colnames(Z) <- allvars
  
  n <- ncol(Y)
  keep.tax <- which(rowSums(Y > 0) / n >= prev_cut)
  Y <- Y[keep.tax, ]
  m <- nrow(Y)
  
  ## some samples may have zero total counts after screening taxa
  if (any(colSums(Y) == 0)) {
    ind <- which(colSums(Y) > 0)
    Y <- Y[, ind]
    Z <- as.data.frame(Z[ind, ])
    names(Z) <- allvars
    keep.sam <- keep.sam[ind]
    n <- ncol(Y)
  }
  
  ## scaling numerical variables
  ind <- sapply(seq_len(ncol(Z)), function(i) is.numeric(Z[, i]))
  Z[, ind] <- scale(Z[, ind])
  p <- ncol(Z) + 1
  
  ## handling zeros
  if (any(Y == 0)) {
    N <- colSums(Y)
    if (adaptive) {
      logN <- log(N)
      tmp <- lm(as.formula(paste0("logN", formula)), Z)
      corr.pval <- coef(summary(tmp))[-1, "Pr(>|t|)"]
      if (any(corr.pval <= corr_cut)) {
        # cat("Imputation approach is used.\n")
        imputation <- TRUE
      } else {
        # cat("Pseudo-count approach is used.\n")
        imputation <- FALSE
      }
    }
    if (imputation) {
      N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
      N.mat[Y > 0] <- 0
      tmp <- N[max.col(N.mat)]
      Y <- Y + N.mat / tmp
    } else {
      Y <- Y + pseudo_cnt
    }
  }
  
  ## CLR transformation
  logY <- log2(Y)
  W <- t(logY) - colMeans(logY)
  
  ## quantile regression
  options(warn = -1)
  formula_use <- as.formula(paste0("w_tmp", formula))
  formula_train <- as.formula(paste0("w_train", formula))
  alpha_vec <- numeric(m)
  sd_alpha <- numeric(m)
  cv_fold <- origami::folds_vfold(n = n, V = n_cv)
  for (iter_m in seq_len(m)) {
    ## sample index
    w_tmp <- W[, iter_m]
    resid_para <- numeric(n_tau)
    resid_cv <- numeric(n_cv)
    for (iter_tau in seq_len(n_tau)) {
      tau_use <- tau_vec[iter_tau]
      for (iter_fold in seq_len(n_cv)) {
        ## load index
        cv_tmp <- cv_fold[[iter_fold]]
        index_train <- cv_tmp$training_set
        index_cross <- cv_tmp$validation_set
        ## load data
        Z_cross <- data.frame(Z[index_cross, ])
        colnames(Z_cross) <- colnames(Z)
        w_cross <- w_tmp[index_cross]
        w_train <- w_tmp[index_train]
        Z_train <- data.frame(Z[index_train, ])
        colnames(Z_train) <- colnames(Z)
        ## fit model
        fit_tmp <- quantreg::rq(formula_train, tau = tau_use, data = Z_train)
        resid_tmp <- w_cross - predict(fit_tmp, data = Z_cross)
        resid_cv[iter_fold] <- sum(resid_tmp * resid_tmp)
      }
      resid_para[iter_tau] <- mean(resid_cv)
    }
    ## fit model
    tau <- tau_vec[which.min(resid_para)]
    fit <- quantreg::rq(formula_use, tau = tau, data = Z)
    alpha_vec[iter_m] <- coef(fit)["u1"]
    summary_fit <- summary(fit, se = "boot", R = 500)
    sd_alpha[iter_m] <- coef(summary_fit)["u1", 2]
  }
  options(warn = 0)
  
  #### mode correction
  # bias <- median(alpha_vec)
  bias <- modeest::mlv(sqrt(n) * alpha_vec,
                       method = "meanshift", kernel = "gaussian"
  ) / sqrt(n)
  alpha_correct <- alpha_vec - bias
  
  #### test
  if (test_method == "t") {
    ## t-test
    stat <- alpha_correct / sd_alpha
    p_value <- 2 * pt(-abs(stat), df = n - p)
    p_adj <- p.adjust(p_value, method = "BH")
    index_select <- which(p_adj < 0.05)
  } else if (test_method == "chi") {
    ## chi square test
    stat <- alpha_correct * alpha_correct / (sd_alpha * sd_alpha)
    p_value <- 1 - pchisq(stat, df = 1)
    p_adj <- p.adjust(p_value, method = "BH")
    index_select <- which(p_adj < 0.05)
  }
  
  #### return results
  return(list(
    alpha_correct = alpha_correct, sd_alpha = sd_alpha, p_value = p_value,
    index_select = index_select
  ))
}

####################### winsorization #######################
winsor.fun <- function(Y, quan) {
  N <- colSums(Y)
  P <- t(t(Y) / N)
  cut <- apply(P, 1, quantile, quan)
  Cut <- matrix(rep(cut, ncol(Y)), nrow(Y))
  ind <- P > Cut
  P[ind] <- Cut[ind]
  Y <- round(t(t(P) * N))
  return(Y)
}

# ####################### winsorization  LinDA #######################
# linda_winsor <- function(otu.tab, meta, formula, type = 'count',
#                          adaptive = TRUE, imputation = FALSE, pseudo.cnt = 0.5, corr.cut = 0.1,
#                          p.adj.method = 'BH', alpha = 0.05,
#                          prev.cut = 0, lib.cut = 1, winsor.quan.vec = NULL, n.cores = 1) {
#   if(any(is.na(otu.tab))) {
#     stop('The OTU table contains NAs! Please remove!\n')
#   }
#   allvars <- all.vars(as.formula(formula))
#   Z_all <- as.data.frame(meta[, allvars])
#   
#   ## preprocessing
#   keep.sam <- which(colSums(otu.tab) >= lib.cut & rowSums(is.na(Z_all)) == 0)
#   Y_all <- otu.tab[, keep.sam]
#   Z_all <- as.data.frame(Z_all[keep.sam, ])
#   names(Z_all) <- allvars
#   
#   n <- ncol(Y_all)
#   keep.tax <- which(rowSums(Y_all > 0) / n >= prev.cut)
#   Y_all <- Y_all[keep.tax, ]
#   m <- nrow(Y_all)
#   
#   ## some samples may have zero total counts after screening taxa
#   if(any(colSums(Y_all) == 0)) {
#     ind <- which(colSums(Y_all) > 0)
#     Y_all <- Y_all[, ind]
#     Z_all <- as.data.frame(Z_all[ind, ])
#     names(Z_all) <- allvars
#     keep.sam <- keep.sam[ind]
#     n <- ncol(Y_all)
#   }
#   
#   ## scaling numerical variables
#   ind <- sapply(1 : ncol(Z_all), function(i) is.numeric(Z_all[, i]))
#   Z_all[, ind] <- scale(Z_all[, ind])
#   
#   ## winsorization
#   if (is.null(winsor.quan.vec)) {
#     winsor.quan.vec <- seq(0.9, 1, by = 0.01)
#     n.winsor <- length(winsor.quan.vec)
#   }
#   
#   fit_list <- list()
#   Y_list <- list()
#   Z_list <- list()
#   W_list <- list()
#   resid_vec <- numeric(n.winsor)
#   for (iter.winsor in seq_len(n.winsor)) {
#     ## winsorization
#     winsor.quan <- winsor.quan.vec[iter.winsor]
#     dta_winsor <- winsor_fun(Y_all, winsor.quan)
#     Y <- dta_winsor$Y
#     n_keep <- dta_winsor$n_keep
#     
#     ## select samples
#     N_tmp <- colSums(Y)
#     keep_sam <- which(N_tmp >= 1)
#     Y <- Y[, keep_sam]
#     allvars <- colnames(Z_all)
#     Z <- as.data.frame(Z[keep_sam, ])
#     colnames(Z) <- allvars
#     n_keep <- n_keep - length(which(N_tmp == 0)) * m
#     
#     ## 
#     if(grepl('\\(', formula)) {
#       random.effect <- TRUE
#     } else {
#       random.effect <- FALSE
#     }
#     
#     if(is.null(rownames(otu.tab))) {
#       taxa.name <- (1 : nrow(otu.tab))[keep.tax]
#     } else {
#       taxa.name <- rownames(otu.tab)[keep.tax]
#     }
#     if(is.null(rownames(meta))) {
#       samp.name <- (1 : nrow(meta))[keep.sam]
#     } else {
#       samp.name <- rownames(meta)[keep.sam]
#     }
#     
#     ## handling zeros
#     if(type == 'count') {
#       if(any(Y == 0)) {
#         N <- colSums(Y)
#         if(adaptive) {
#           logN <- log(N)
#           if(random.effect) {
#             tmp <- lmer(as.formula(paste0('logN', formula)), Z)
#           } else {
#             tmp <- lm(as.formula(paste0('logN', formula)), Z)
#           }
#           corr.pval <- coef(summary(tmp))[-1, "Pr(>|t|)"]
#           if(any(corr.pval <= corr.cut)) {
#             imputation <- TRUE
#           } else {
#             imputation <- FALSE
#           }
#         }
#         if(imputation) {
#           N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
#           N.mat[Y > 0] <- 0
#           tmp <- N[max.col(N.mat)]
#           Y <- Y + N.mat / tmp
#         } else {
#           Y <- Y + pseudo.cnt
#         }
#       }
#     }
#     
#     if(type == 'proportion') {
#       if(any(Y == 0)) {
#         ## Half minimum approach
#         Y <- t(apply(Y, 1, function (x) {
#           x[x == 0] <- 0.5 * min(x[x != 0])
#           return(x)
#         }))
#       }
#     }
#     
#     ## CLR transformation
#     logY <- log2(Y)
#     W <- t(logY) - colMeans(logY)
#     Y_list[[iter.winsor]] <- Y
#     Z_list[[iter.winsor]] <- Z
#     W_list[[iter.winsor]] <- W
#     
#     
#     ## linear regression
#     oldw <- getOption('warn')
#     options(warn = -1)
#     if (!random.effect) {
#       suppressMessages(fit <- lm(as.formula(paste0('W', formula)), Z))
#       fit_list[[iter.winsor]] <- fit
#       resid_all <- t(resid(fit))
#       resid_vec[iter.winsor] <- sum(resid_all * resid_all) / (n_keep * (n_keep - 1))
#     }
#   }
#   
#   #### select the best quantile
#   iter_select <- which.min(resid_vec)
#   winsor.quan <- winsor.quan.vec[iter_select]
#   print(winsor.quan)
#   W <- W_list[[iter_select]]
#   Y <- Y_list[[iter_select]]
#   Z <- Z_list[[iter_select]]
#   fit <- fit_list[[iter_select]]
#   
#   if(!random.effect) {
#     res <- do.call(rbind, coef(summary(fit)))
#     d <- ncol(model.matrix(fit))
#     df <- rep(n - d, m)
#     tmp <- vcov(fit)
#     res.cov <- lapply(seq_len(m), function(i) {
#       tmp[((i-1)*d+1) : (i*d), ((i-1)*d+1) : (i*d)]
#     })
#     res.cov <- do.call(rbind, res.cov)
#     rownames(res.cov) <- rownames(res)
#     colnames(res.cov) <- rownames(res)[1 : d]
#   } else {
#     fun <- function(i) {
#       w <- W[, i]
#       fit <- lmer(as.formula(paste0('w', formula)), Z)
#       list(coef(summary(fit)), vcov(fit))
#     }
#     if(n.cores > 1) {
#       tmp <- mclapply(c(1 : m), function(i) fun(i), mc.cores = n.cores)
#     } else {
#       suppressMessages(tmp <- lapply(seq_len(m), function(i) {return(fun(i))}))
#     }
#     res <- do.call(rbind, lapply(tmp, `[[`, 1))
#     res.cov <- do.call(rbind, lapply(tmp, `[[`, 2))
#   }
#   options(warn = oldw)
#   
#   res.intc <- res[which(rownames(res) == '(Intercept)'), ]
#   rownames(res.intc) <- NULL
#   options(warn = -1)
#   suppressMessages(tmp <- modeest::mlv(sqrt(n) * res.intc[, 1],
#                                        method = 'meanshift', kernel = 'gaussian') / sqrt(n))
#   options(warn = oldw)
#   baseMean <- 2 ^ (res.intc[, 1] - tmp)
#   baseMean <- baseMean / sum(baseMean) * 1e6
#   
#   output.fun <- function(x) {
#     res.voi <- res[which(rownames(res) == x), ]
#     rownames(res.voi) <- NULL
#     
#     if(random.effect) {
#       df <- res.voi[, 3]
#     }
#     
#     log2FoldChange <- res.voi[, 1]
#     lfcSE <- res.voi[, 2]
#     oldw <- getOption('warn')
#     options(warn = -1)
#     suppressMessages(bias <- modeest::mlv(sqrt(n) * log2FoldChange,
#                                           method = 'meanshift', kernel = 'gaussian') / sqrt(n))
#     options(warn = oldw)
#     log2FoldChange <- log2FoldChange - bias
#     stat <- log2FoldChange / lfcSE
#     
#     pvalue <- 2 * pt(-abs(stat), df)
#     padj <- p.adjust(pvalue, method = p.adj.method)
#     reject <- padj <= alpha
#     output <- cbind.data.frame(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, reject, df)
#     rownames(output) <- taxa.name
#     return(list(bias = bias, output = output))
#   }
#   
#   cov.fun <- function(x) {
#     tmp <- (1 : ncol(res.cov))[-c(1, which(colnames(res.cov) == x))]
#     covariance <- as.data.frame(as.matrix(res.cov[which(rownames(res.cov) == x), tmp]))
#     rownames(covariance) <- taxa.name
#     colnames(covariance) <- colnames(res.cov)[tmp]
#     return(covariance)
#   }
#   
#   variables <- unique(rownames(res))[-1]
#   variables.n <- length(variables)
#   bias <- rep(NA, variables.n)
#   output <- list()
#   if(variables.n == 1) {
#     covariance <- NULL
#   } else {
#     covariance <- list()
#   }
#   for(i in 1 : variables.n) {
#     tmp <- output.fun(variables[i])
#     output[[i]] <- tmp[[2]]
#     bias[i] <- tmp[[1]]
#     if(variables.n > 1) {
#       covariance[[i]] <- cov.fun(variables[i])
#     }
#   }
#   names(output) <- variables
#   if(variables.n > 1) {
#     names(covariance) <- variables
#   }
#   
#   rownames(Y) <- taxa.name
#   colnames(Y) <- samp.name
#   rownames(Z) <- samp.name
#   return(list(variables = variables, bias = bias, output = output, covariance = covariance, otu.tab.use = Y, meta.use = Z))
# }