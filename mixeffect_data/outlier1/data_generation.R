set.seed(1)
library(parallel)
## load parameters
para0 <- readRDS("mixeffect_data/outlier1/datasets/log.normal.para.rds")
beta0 <- para0$beta0
sigma2 <- para0$sigma2
# parameter use
sample.size.vec <- 200
m <- 500
n_sim <- 100
# define settings
sig.density.vec <- c(0.05, 0.2)
sig.strength.vec <- seq(1.05, 2, length.out = 6)
s1 <- 1
s2 <- 2
s3 <- 6
sample.size <- rep(sample.size.vec, each = s2 * s3)
sig.density <- rep(rep(sig.density.vec, each = s3), s1)
sig.strength <- rep(sig.strength.vec, s1 * s2)
setting <- cbind(sample.size, sig.density, sig.strength)
colnames(setting) <- c("sample.size", "sig.density", "sig.strength")
n_setting <- nrow(setting)

#### generate datasets
for (iter_para in seq_len(n_setting)) {
  print(iter_para)
  para <- setting[iter_para, ]
  n <- as.numeric(para[1])
  gamma <- as.numeric(para[2])
  mu_use <- as.numeric(para[3])
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

    if (n == 50) {
      n.id <- 25
      u <- c(rep(0, 24), rep(1, 26))
      id <- rep(1:n.id, each = 2)
    } else if (n == 200) {
      n.id <- 50
      u <- c(rep(0, 100), rep(1, 100))
      id <- rep(1:n.id, each = 4)
    }
    ## t random effect
    r <- matrix(rt(m * n.id, df = 4), nrow = m)[, id]
    Z <- cbind(u, id)
    beta <- alpha
    tmp <- beta %*% t(Z[, 1]) + beta0
    logX <- tmp + matrix(rt(m * n, df = 3), nrow = m) + r
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
    "mixeffect_data/outlier1/datasets/n", n,
    "gamma", gamma, "mu", mu_use, ".rds"
  ))
}
