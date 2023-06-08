####################### robust regression #######################
rlm_fun <- function(Y, Z, formula, adaptive = TRUE, imputation = FALSE,
                    pseudo_cnt = 0.5, corr_cut = 0.1, prev_cut = 0, lib_cut = 1,
                    res_method = "psi.huber", test_method = "t",
                    adj_method = "BH", alpha = 0.05) {
  #### check input
  if (any(is.na(Y))) {
    stop("The OTU table contains NAs! Please remove!\n")
  }
  allvars <- colnames(Z)
  
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
  
  ## robust linear regression
  options(warn = -1)
  formula_use <- as.formula(paste0("w_tmp", formula))
  formula_train <- as.formula(paste0("w_train", formula))
  alpha_vec <- numeric(m)
  sd_alpha <- numeric(m)
  
  ## Huber method
  if (res_method == "psi.huber") {
    hyper_para_vec <- exp(seq(log(0.2), log(5), length = 20))
    cv_fold <- origami::folds_vfold(n = n, V = 5)
    for (iter_m in seq_len(m)) {
      ## sample index
      w_tmp <- W[, iter_m]
      resid_para <- numeric(20)
      resid_cv <- numeric(5)
      for (iter_cross in seq_len(20)) {
        hyper_para_use <- hyper_para_vec[iter_cross]
        for (iter_fold in seq_len(5)) {
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
          fit_tmp <-  MASS::rlm(formula_train, data = Z_train, psi = res_method, maxit = 50,
                                k = hyper_para_use)
          resid_cv[iter_fold] <- sum((w_cross - predict(fit_tmp, data = Z_cross))^2)
        } 
        resid_para[iter_cross] <- mean(resid_cv)
      }
      ## fit model
      hyper_para <- hyper_para_vec[which.min(resid_para)]
      fit <- MASS::rlm(formula_use, data = Z, psi = res_method, maxit = 50, k = hyper_para)
      alpha_vec[iter_m] <- coef(fit)["u1"]
      summary_fit <- summary(fit)
      sd_alpha[iter_m] <- coef(summary_fit)["u1", 2]
    }
  } else if (res_method == "psi.bisquare") {
    hyper_para_vec <- exp(seq(log(1), log(10), length = 20))
    cv_fold <- origami::folds_vfold(n = n, V = 5)
    for (iter_m in seq_len(m)) {
      ## sample index
      w_tmp <- W[, iter_m]
      resid_para <- numeric(20)
      resid_cv <- numeric(5)
      for (iter_cross in seq_len(20)) {
        hyper_para_use <- hyper_para_vec[iter_cross]
        for (iter_fold in seq_len(5)) {
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
          fit_tmp <-  MASS::rlm(formula_train, data = Z_train, psi = res_method, maxit = 50,
                                c = hyper_para_use)
          resid_cv[iter_fold] <- sum((w_cross - predict(fit_tmp, data = Z_cross))^2)
        } 
        resid_para[iter_cross] <- mean(resid_cv)
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
qr_fun <- function(Y, Z, formula, tau = NULL, adaptive = TRUE, imputation = FALSE,
                   pseudo_cnt = 0.5, corr_cut = 0.1, prev_cut = 0, lib_cut = 1,
                   test_method = "t", adj_method = "BH", alpha = 0.05) {
  #### check input
  if (any(is.na(Y))) {
    stop("The OTU table contains NAs! Please remove!\n")
  }
  allvars <- colnames(Z)
  if (is.null(tau)) {
    tau <- 0.5
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
  alpha_vec <- numeric(m)
  sd_alpha <- numeric(m)
  for (iter_m in seq_len(m)) {
    w_tmp <- W[, iter_m]
    fit <- quantreg::rq(formula_use, tau = tau, data = Z)
    alpha_vec[iter_m] <- coef(fit)["u1"]
    summary_fit <- summary(fit, se = "boot", R = 500)
    sd_alpha[iter_m] <- coef(summary_fit)["u1", 2]
  }
  options(warn = 0)
  
  #### mode correction
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
