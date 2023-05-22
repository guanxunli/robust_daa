set.seed(2023)
library(parallel)
## load parameters
para0 <- readRDS("loglinear_data/outlier2/datasets/log.normal.para.rds")
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

#### without confounder
for (iter_para in seq_len(nrow(setting))) {
  print(iter_para)
  para <- setting[iter_para, ]
  n <- para[1]
  gamma <- para[2]
  mu_use <- para[3]
  dta_list <- mclapply(seq_len(n_sim), function(iter_simu) {
    dta <- list()
    ## generate basic expression
    X0 <- matrix(exp(rnorm(m * n, beta0, sqrt(sigma2))), nrow = m)
    pi0 <- t(t(X0) / colSums(X0))
    pi0.ave <- rowMeans(pi0)
    if (any(pi0.ave == 0)) {
      ind <- which(pi0.ave == 0)
      pi0.ave[ind] <- min(pi0.ave[-ind]) / 10
      pi0.ave <- pi0.ave / sum(pi0.ave)
    }
    tmp <- (pi0.ave > 0.005)
    mu <- 2 * mu_use * (n <= 50) + mu_use * (n > 50)
    mu.1 <- log(mu * tmp + mu * (0.005 / pi0.ave)^(1 / 3) * (1 - tmp))

    ## index true
    index_true <- rbinom(m, 1, gamma)
    index_alter <- which(index_true == 1)
    alpha <- mu.1 * index_true

    ## generative confounders
    u <- rbinom(n, 1, 0.5)
    Z <- cbind(u)
    beta <- alpha

    ## generate X
    tmp <- beta %*% t(Z) + beta0
    error_mat <- matrix(rlnorm(m * n, meanlog = log(1.2), sdlog = log(2)), nrow = m)
    error_mat <- error_mat - mean(error_mat)

    logX <- tmp + error_mat
    pi <- apply(logX, 2, function(logx) {
      max_logx <- max(logx)
      x <- exp(logx - max_logx)
      return(x / sum(x))
    })
    # generate Y
    N <- rnbinom(n, size = 5.3, mu = 7645)
    Y <- sapply(1:n, function(s) rmultinom(1, N[s], pi[, s]))

    ## save results
    Z <- as.data.frame(Z)
    Z$u <- factor(Z$u)
    Y <- as.data.frame(Y)
    colnames(Y) <- rownames(Z) <- paste0("sample", 1:n)
    rownames(Y) <- paste0("taxon", 1:m)
    dta$Y <- Y
    dta$Z <- Z
    dta$index_alter <- index_alter
    return(dta)
  }, mc.cores = 50)
  # save datasets
  saveRDS(dta_list, paste0(
    "loglinear_data/outlier2/datasets/nocon_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))
}

#### with confounder
for (iter_para in seq_len(nrow(setting))) {
  print(iter_para)
  para <- setting[iter_para, ]
  n <- para[1]
  gamma <- para[2]
  mu_use <- para[3]
  dta_list <- mclapply(seq_len(n_sim), function(iter_simu) {
    dta <- list()
    ## generate basic expression
    X0 <- matrix(exp(rnorm(m * n, beta0, sqrt(sigma2))), nrow = m)
    pi0 <- t(t(X0) / colSums(X0))
    pi0.ave <- rowMeans(pi0)
    if (any(pi0.ave == 0)) {
      ind <- which(pi0.ave == 0)
      pi0.ave[ind] <- min(pi0.ave[-ind]) / 10
      pi0.ave <- pi0.ave / sum(pi0.ave)
    }
    tmp <- (pi0.ave > 0.005)
    mu <- 2 * mu_use * (n <= 50) + mu_use * (n > 50)
    mu.1 <- log(mu * tmp + mu * (0.005 / pi0.ave)^(1 / 3) * (1 - tmp))

    ## index true
    index_true <- rbinom(m, 1, gamma)
    index_alter <- which(index_true == 1)
    alpha <- mu.1 * index_true

    ## generative confounders
    # c mat
    z1 <- rbinom(n, 1, 0.5)
    z1[which(z1 == 0)] <- -1
    z2 <- rnorm(n, 0, 1)
    beta1 <- rnorm(m, 1, 1)
    beta2 <- rnorm(m, 2, 1)
    # u mat
    u <- rbinom(n, 1, 1 / (1 + exp(-0.5 * z1 - 0.5 * z2)))
    while (length(unique(u)) == 1) {
      u <- rbinom(n, 1, 1 / (1 + exp(-z1 - z2)))
    }
    Z <- cbind(u, z1, z2)
    beta <- cbind(alpha, beta1, beta2)

    ## generate X
    tmp <- beta %*% t(Z) + beta0
    # error mat
    error_mat <- matrix(rlnorm(m * n, meanlog = log(1.2), sdlog = log(2)), nrow = m)
    error_mat <- error_mat - mean(error_mat)
    # log X
    logX <- tmp + error_mat
    pi <- apply(logX, 2, function(logx) {
      max_logx <- max(logx)
      x <- exp(logx - max_logx)
      return(x / sum(x))
    })
    # generate Y
    N <- rnbinom(n, size = 5.3, mu = 7645)
    Y <- sapply(1:n, function(s) rmultinom(1, N[s], pi[, s]))

    ## save results
    Z <- as.data.frame(Z)
    Z$u <- as.factor(Z$u)
    Z$z1 <- as.factor(Z$z1)
    Y <- as.data.frame(Y)
    colnames(Y) <- rownames(Z) <- paste0("sample", 1:n)
    rownames(Y) <- paste0("taxon", 1:m)
    dta$Y <- Y
    dta$Z <- Z
    dta$index_alter <- index_alter
    return(dta)
  }, mc.cores = 50)
  # save datasets
  saveRDS(dta_list, paste0(
    "loglinear_data/outlier2/datasets/con_n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))
}
