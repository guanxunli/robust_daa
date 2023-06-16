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

####################### Cauchy method #######################
CombinePValues <- function(pvalues, weights = NULL) {
  if (!is.matrix(pvalues)) {
    pvalues <- as.matrix(pvalues)
  }
  ## to avoid extremely values
  pvalues[which(pvalues == 0)] <- 5.55e-17
  pvalues[which((1 - pvalues) < 1e-3)] <- 0.99
  
  num_pval <- ncol(pvalues)
  num_gene <- nrow(pvalues)
  if (is.null(weights)) {
    weights <- matrix(rep(1.0 / num_pval, num_pval * num_gene), ncol = num_pval)
  } # end fi
  if ((nrow(weights) != num_gene) || (ncol(weights) != num_pval)) {
    stop("the dimensions of weights does not match that of combined pvalues")
  } # end fi
  
  Cstat <- tan((0.5 - pvalues) * pi)
  
  wCstat <- weights * Cstat
  Cbar <- apply(wCstat, 1, sum)
  # combined_pval <- 1.0/2.0 - atan(Cbar)/pi
  combined_pval <- 1.0 - pcauchy(Cbar)
  combined_pval[which(combined_pval <= 0)] <- 5.55e-17
  return(combined_pval)
}

####################### LinDA winsorization #######################
linda_winsor <- function(Y, Z, formula, winsor.quan.vec = NULL) {
  ## initialization
  m <- nrow(Y)
  if (is.null(winsor.quan.vec)) {
    winsor.quan.vec <- seq(0.9, 1, by = 0.01)
  }
  n_winsor <- length(winsor.quan.vec)
  p_value_mat <- matrix(0, nrow = m, ncol = n_winsor)
  index_rm <- NULL
  ## different winsorization
  for (iter.winsor in seq_len(n_winsor)) {
    ## initialization
    winsor.quan <- winsor.quan.vec[iter.winsor]
    Y_use <- Y
    Z_use <- Z
    ## remove samples for Y
    Y_tmp <- winsor.fun(Y, winsor.quan)
    N_tmp <- colSums(Y_tmp)
    keep_sam <- which(N_tmp >= 1)
    Y_use <- Y_use[, keep_sam]
    ## remove samples for Z
    allvars <- colnames(Z)
    Z_use <- as.data.frame(Z_use[keep_sam, ])
    colnames(Z_use) <- allvars
    ## linda method
    res <- try(LinDA::linda(
      otu.tab = Y_use, meta = Z_use, formula = formula,
      winsor.quan = winsor.quan
    ), silent = TRUE)
    if (class(res) == "try-error") {
      index_rm <- c(index_rm, iter.winsor)
    } else {
      p_value_mat[, iter.winsor] <- res$output[[1]]$pvalue
    }
  }
  if (length(index_rm) > 0) {
    p_value_mat <- p_value_mat[, -index_rm]
  }
  pcol <- ncol(p_value_mat)
  weights <- matrix(rep(1.0 / pcol, pcol * m), ncol = pcol)
  combine_pvalue <- CombinePValues(p_value_mat, weights)
  p_adj <- p.adjust(combine_pvalue, "BH")
  index_select <- which(p_adj < 0.05)
  return(list(
    p_value_mat = p_value_mat, combine_pvalue = combine_pvalue,
    index_select = index_select
  ))
}

####################### Huber method #######################
library(lme4)
library(robustlmm)
rlm_fun <- function(Y, Z, formula, adaptive = TRUE, imputation = FALSE,
                    pseudo_cnt = 0.5, corr_cut = 0.1, prev_cut = 0, lib_cut = 1,
                    hyper_para_vec = NULL,
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
  if (is.null(hyper_para_vec)) {
    hyper_para_vec <- exp(seq(log(1.345), log(4), length = 5))
  }
  n_hyper <- length(hyper_para_vec)
  pvalue_mat <- matrix(NA, nrow = m, ncol = n_hyper)
  alpha_mat <- matrix(NA, nrow = m, ncol = n_hyper)
  sd_mat <- matrix(NA, nrow = m, ncol = n_hyper)
  df_mat <- matrix(NA, nrow = m, ncol = n_hyper)
  for (iter_m in seq_len(m)) {
    ## fit robust rlm
    w_tmp <- W[, iter_m]
    rfit <- suppressMessages(rlmer(formula_use, data = Z))
    ## parameter
    summary_rfit <- summary(rfit)
    alpha_mat[iter_m, 1] <- coef(summary_rfit)[2, 1]
    sd_mat[iter_m, 1] <- coef(summary_rfit)[2, 2]
    ## degree of freedom
    rfit_df <- rfit
    class(rfit_df) <- "lmerMod" # to use get_Lb_ddf
    df_KR <- pbkrtest::get_Lb_ddf(rfit_df, lme4::fixef(rfit_df))
    df_mat[iter_m, 1] <- df_KR
    ## different hyperparameter
    for (iter_hyper in c(2:n_hyper)) {
      hyper_para <- hyper_para_vec[iter_hyper]
      rfit_hyper <- update(rfit, rho.sigma.e = psi2propII(smoothPsi, k = hyper_para),
                           rho.sigma.b = psi2propII(smoothPsi, k = hyper_para))
      ## parameter
      summary_rfit <- summary(rfit_hyper)
      alpha_mat[iter_m, iter_hyper] <- coef(summary_rfit)[2, 1]
      sd_mat[iter_m, iter_hyper] <- coef(summary_rfit)[2, 2]
      ## degree of freedom
      rfit_df <- rfit_hyper
      class(rfit_df) <- "lmerMod" # to use get_Lb_ddf
      df_KR <- pbkrtest::get_Lb_ddf(rfit_df, lme4::fixef(rfit_df))
      df_mat[iter_m, iter_hyper] <- df_KR
    }
  }
  options(warn = 0)
  
  for (iter_hyper in seq_len(n_hyper)) {
    alpha_vec <- alpha_mat[, iter_hyper]
    sd_alpha <- sd_mat[, iter_hyper]
    df_vec <- df_mat[, iter_hyper]
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
    } else if (test_method == "chi") {
      ## chi square test
      stat <- alpha_correct * alpha_correct / (sd_alpha * sd_alpha)
      p_value <- 1 - pchisq(stat, df = 1)
    }
    pvalue_mat[, iter_hyper] <- p_value
  }
  weights <- matrix(rep(1.0 / n_hyper, n_hyper * m), ncol = n_hyper)
  combine_pvalue <- CombinePValues(pvalue_mat, weights)
  p_adj <- p.adjust(combine_pvalue, "BH")
  index_select <- which(p_adj < 0.05)
  #### return results
  return(list(
    p_value_mat = pvalue_mat, combine_pvalue = combine_pvalue,
    index_select = index_select
  ))
}