set.seed(1)
library(parallel)
#### load dataset
dta <- readRDS("real_data/datasets/CDI_IBD_RA_SMOKE.rds")
preprocess_fun <- function(otu_tab, meta, prev_cut = 0.1, lib.cut = 1000) {
  keep_sam <- colSums(otu_tab) >= lib.cut
  Y <- otu_tab[, keep_sam]
  Z <- as.data.frame(meta[keep_sam, ])
  names(Z) <- names(meta)
  rownames(Z) <- rownames(meta)[keep_sam]
  keep_tax <- rowSums(Y > 0) / ncol(Y) >= prev_cut
  Y <- Y[keep_tax, ]
  return(list(Y = Y, Z = Z, keep_sam = keep_sam, keep_tax = keep_tax))
}
thres_vec <- c(0.01, seq(0.05, 0.25, 0.05))
n_thres <- length(thres_vec)

#### Run different methods
formulas <- c('Disease', 'Disease+Antibiotic', 'Disease')
n_permu <- 1e3
for (iter_dta in seq_len(3)) {
  ## load data set
  otu_tab <- dta[[2 * (iter_dta - 1) + 1]]
  meta <- dta[[2 * iter_dta]]
  dta_process <- preprocess_fun(otu_tab, meta)
  otu_tab <- dta_process$Y
  meta <- dta_process$Z
  formula <- formulas[iter_dta]
  allvars <- colnames(meta)
  out_res_list <- mclapply(seq_len(n_permu), function (iter_permu) {
    meta_permu <- meta
    meta_permu[, allvars] <- sample(meta[, allvars])
    out_res <- list()
    ## LinDA method
    linda_res <- LinDA::linda(otu_tab, meta_permu, paste("~", formula))
    out_res[["linda"]] <- which(linda_res$output[[1]]$reject == TRUE)
    
    ## LinDA97
    ## remove samples for Y
    otu_tab_tmp <- robustDAA::winsor.fun(otu_tab, 0.97)
    N_tmp <- colSums(otu_tab_tmp)
    keep_sam <- which(N_tmp >= 1)
    otu_tab_use <- otu_tab[, keep_sam]
    ## remove samples for Z
    meta_use <- as.data.frame(meta_permu[keep_sam, ])
    colnames(meta_use) <- allvars
    ## linda method
    linda97_res <- LinDA::linda(
      otu.tab = otu_tab_use, meta = meta_use, formula = paste("~", formula),
      winsor.quan = 0.97
    )
    out_res[["linda97"]] <- which(linda97_res$output[[1]]$reject == TRUE)
    
    ## LinDA90
    ## remove samples for Y
    otu_tab_tmp <- robustDAA::winsor.fun(otu_tab, 0.90)
    N_tmp <- colSums(otu_tab_tmp)
    keep_sam <- which(N_tmp >= 1)
    otu_tab_use <- otu_tab[, keep_sam]
    ## remove samples for Z
    meta_use <- as.data.frame(meta_permu[keep_sam, ])
    colnames(meta_use) <- allvars
    ## linda method
    linda90_res <- LinDA::linda(
      otu.tab = otu_tab_use, meta = meta_use, formula = paste("~", formula),
      winsor.quan = 0.90
    )
    out_res[["linda90"]] <- which(linda90_res$output[[1]]$reject == TRUE)
    
    ## Huber method
    huber_res <- robustDAA::rlm_fun(
      otu_tab = otu_tab, meta = meta_permu, formula = paste("~", formula),
      reg_method = "psi.huber"
    )
    out_res[["huber"]] <- huber_res$index_select
    
    ## Bisquare method
    bisquareres <- robustDAA::rlm_fun(
      otu_tab = otu_tab, meta = meta_permu, formula = paste("~", formula),
      reg_method = "psi.bisquare"
    )
    out_res[["bisquare"]] <- bisquareres$index_select
    
    ## QR method
    qrres <- robustDAA::qr_fun(otu_tab = otu_tab, meta = meta_permu, paste("~", formula))
    out_res[["qr"]] <- qrres$index_select
    
    ## return results
    return(out_res)
  }, mc.cores = 50)
  saveRDS(out_res_list, paste0("real_data/results/data", iter_dta, "res.rds"))
}