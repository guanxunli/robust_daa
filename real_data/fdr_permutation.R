library(parallel)
#### load dataset
dta <- readRDS("real_data/datasets/CDI_IBD_RA_SMOKE.rds")
preprocess_fun <- function(otu_tab, meta, prev_cut = 0.1, lib_cut = 1000) {
  keep_sam <- colSums(otu_tab) >= lib_cut
  Y <- otu_tab[, keep_sam]
  Z <- as.data.frame(meta[keep_sam, ])
  names(Z) <- names(meta)
  rownames(Z) <- rownames(meta)[keep_sam]
  keep_tax <- rowSums(Y > 0) / ncol(Y) >= prev_cut
  Y <- Y[keep_tax, ]
  return(list(Y = Y, Z = Z, keep_sam = keep_sam, keep_tax = keep_tax))
}

#### Run different methods
formulas <- c("Disease", "Disease+Antibiotic", "Disease")
n_permu <- 1e3
for (iter_dta in seq_len(3)) {
  set.seed(1)
  ## load data set
  otu_tab <- dta[[2 * (iter_dta - 1) + 1]]
  meta <- dta[[2 * iter_dta]]
  dta_process <- preprocess_fun(otu_tab, meta)
  otu_tab <- dta_process$Y
  meta <- dta_process$Z
  formula <- formulas[iter_dta]
  allvars <- colnames(meta)
  n_sample <- nrow(meta)
  out_res_list <- mclapply(seq_len(n_permu), function(iter_permu) {
    meta_permu <- as.data.frame(meta[sample(seq_len(n_sample)), ])
    colnames(meta_permu) <- allvars
    out_res <- list()
    ## LinDA method
    linda_res <- LinDA::linda(otu_tab, meta_permu, paste("~", formula))
    padj_linda <- linda_res$output[[1]]$padj
    out_res[["linda"]] <- padj_linda

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
    padj_linda97 <- linda97_res$output[[1]]$padj
    out_res[["linda97"]] <- padj_linda97

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
    padj_linda90 <- linda90_res$output[[1]]$padj
    out_res[["linda90"]] <- padj_linda90

    ## Huber method
    huber_res <- robustDAA::rlm_fun(
      otu_tab = otu_tab, meta = meta_permu, formula = paste("~", formula),
      reg_method = "psi.huber"
    )
    pvalue_huber <- huber_res$combine_pvalue
    padj_huber <- p.adjust(pvalue_huber, method = "BH")
    out_res[["huber"]] <- padj_huber

    ## Bisquare method
    bisquareres <- robustDAA::rlm_fun(
      otu_tab = otu_tab, meta = meta_permu, formula = paste("~", formula),
      reg_method = "psi.bisquare"
    )
    pvalue_bisquare <- bisquareres$combine_pvalue
    padj_bisquare <- p.adjust(pvalue_bisquare, method = "BH")
    out_res[["bisquare"]] <- padj_bisquare

    ## QR method
    qrres <- robustDAA::qr_fun(otu_tab = otu_tab, meta = meta_permu, paste("~", formula))
    pvalue_qr <- qrres$combine_pvalue
    padj_qr <- p.adjust(pvalue_qr, method = "BH")
    out_res[["qr"]] <- padj_qr

    ## return results
    return(out_res)
  }, mc.cores = 50)
  saveRDS(out_res_list, paste0("real_data/results/data", iter_dta, "permu.rds"))
}

#### chekc FDR
method_vec <- c("linda", "linda97", "linda90", "huber", "bisquare", "qr")
n_method <- length(method_vec)
alpha <- 0.1
for (iter_dta in seq_len(3)) {
  cat("Data set", iter_dta, "\n")
  out_res_list <- readRDS(paste0("real_data/results/data", iter_dta, "permu.rds"))
  n_permu <- length(out_res_list)
  out_list <- list()
  for (iter_method in method_vec) {
    out_list[[iter_method]] <- 0
  }
  for (iter_permu in seq_len(n_permu)) {
    out_res <- out_res_list[[iter_permu]]
    for (iter_method in method_vec) {
      padj <- out_res[[iter_method]]
      n_select <- length(which(padj < alpha))
      if (n_select > 0) {
        out_list[[iter_method]] <- out_list[[iter_method]] + 1
      }
    }
  }
  fdr_vec <- as.numeric(unlist(out_list)) / n_permu
  names(fdr_vec) <- method_vec
  print(fdr_vec)
}
