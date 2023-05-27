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

####################### Huber method #######################
rlm_fun <- function(Y, Z, formula, adaptive = TRUE, imputation = FALSE,
                    pseudo_cnt = 0.5, corr_cut = 0.1, prev_cut = 0, lib_cut = 1,
                    test_method = "t", adj_method = "BH", alpha = 0.05) {
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
  
  ## check random effect
  if (grepl("\\(", formula)) {
    random_effect <- TRUE
  } else {
    random_effect <- FALSE
  }
  
  ## handling zeros
  if (any(Y == 0)) {
    N <- colSums(Y)
    if (adaptive) {
      logN <- log(N)
      if (random_effect) {
        tmp <- lmerTest::lmer(as.formula(paste0("logN", formula)), Z)
      } else {
        tmp <- lm(as.formula(paste0("logN", formula)), Z)
      }
      corr_pval <- coef(summary(tmp))[-1, "Pr(>|t|)"]
      if (any(corr_pval <= corr_cut)) {
        # cat("Imputation approach is used.\n")
        imputation <- TRUE
      } else {
        # cat("Pseudo-count approach is used.\n")
        imputation <- FALSE
      }
    }
    if (imputation) {
      Nmat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
      Nmat[Y > 0] <- 0
      tmp <- N[max.col(Nmat)]
      Y <- Y + Nmat / tmp
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
  alpha_vec <- numeric(m)
  sd_alpha <- numeric(m)
  df_vec <- numeric(m)
  if (random_effect) {
    for (iter_m in seq_len(m)) {
      w_tmp <- W[, iter_m]
      rfit <- suppressMessages(robustlmm::rlmer(formula_use, data = Z))
      summary_rfit <- summary(rfit)
      class(rfit) <- "lmerMod" # to use get_Lb_ddf
      df_KR <- pbkrtest::get_Lb_ddf(rfit, lme4::fixef(rfit))
      alpha_vec[iter_m] <- coef(summary_rfit)[2, 1]
      sd_alpha[iter_m] <- coef(summary_rfit)[2, 2]
      df_vec[iter_m] <- df_KR
    }
  } else {
    for (iter_m in seq_len(m)) {
      w_tmp <- W[, iter_m]
      fit <- MASS::rlm(formula_use, data = Z, maxit = 50)
      alpha_vec[iter_m] <- coef(fit)["u1"]
      summary_fit <- summary(fit)
      sd_alpha[iter_m] <- coef(summary_fit)["u1", 2]
    }
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
    p_value <- 2 * pt(-abs(stat), df = df_vec)
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
