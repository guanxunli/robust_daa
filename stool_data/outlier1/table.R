set.seed(1)
outlier <- "outlier1"
################## define parameters ##################
ratio_vec <- paste0("ratio", seq_len(4))
n_ratio <- length(ratio_vec)
method_vec <- c(
  "linda", "linda97", "linda90", "huber", "bisquare", "qr"
)
n_method <- length(method_vec)
## define parameters
n_sam <- 50
n_taxa <- 500
set_df <- data.frame(n_sam = rep(n_sam, length(n_taxa)), n_taxa = rep(n_taxa, each = length(n_sam)))
signa_den <- c(0.05, 0.2)
set_df <- cbind(apply(set_df, 2, rep, length(signa_den)), rep(signa_den, each = nrow(set_df)))
diff_mode <- c("abundant", "mix")
set_df <- cbind(apply(set_df, 2, rep, length(diff_mode)), rep(diff_mode, each = nrow(set_df)))
colnames(set_df) <- c("n_sam", "n_taxa", "signa_den", "mode")
set_df <- as.data.frame(set_df)
nset <- nrow(set_df)
signa_streng <- 1
n_simu <- 1e2
res_list <- list()

################## without confounder ##################
for (iter_set in seq_len(nset)) {
  n_sam <- as.numeric(set_df$n_sam[iter_set])
  n_taxa <- as.numeric(set_df$n_taxa[iter_set])
  signa_den <- as.numeric(set_df$signa_den[iter_set])
  diff_mode <- set_df$mode[iter_set]
  print(c(n_sam, n_taxa, signa_den, diff_mode))
  res_mat <- matrix(NA, nrow = n_method, ncol = 2 * n_ratio)
  #### each ratio
  for (iter_ratio in seq_len(n_ratio)) {
    ratio <- ratio_vec[iter_ratio]
    #### without confounders
    dta_list <- readRDS(paste0(
      "stool_data/", outlier, "/", ratio, "/datasets/noconf_nsam", n_sam, "ntaxa", n_taxa,
      "signal", signa_den, "streng", signa_streng, "mode", diff_mode, ".rds"
    ))
    res_table_noconf <- data.frame(
      "method" = c("LinDA", "LinDA97", "LinDA90", "Huber", "Bi-square", "QR"),
      "confounder" = rep("without", n_method),
      "Power" = numeric(n_method), "Power_sd" = numeric(n_method),
      "FDR" = numeric(n_method), "FDR_sd" = numeric(n_method)
    )
    ## load results
    res_list <- list()
    res_mat_list <- list()
    for (iter_method in method_vec) {
      res_list[[iter_method]] <- readRDS(paste0(
        "stool_data/", outlier, "/", ratio, "/results/", iter_method, "_noconf_nsam", n_sam, "ntaxa", n_taxa,
        "signal", signa_den, "streng", signa_streng, "mode", diff_mode, ".rds"
      ))
      res_mat_list[[iter_method]] <- matrix(NA, nrow = n_simu, ncol = 2)
      colnames(res_mat_list[[iter_method]]) <- c("power", "fdr")
    }
    ## calculate results
    for (iter_simu in seq_len(n_simu)) {
      dta <- dta_list[[iter_simu]]
      index_alter <- which(dta$diff.otu.ind == TRUE)
      for (iter_method in method_vec[seq_len(n_method)]) {
        index_select <- res_list[[iter_method]][[iter_simu]]
        if (length(index_select) == 0) {
          res_mat_list[[iter_method]][iter_simu, ] <- 0
        } else {
          res_mat_list[[iter_method]][iter_simu, 1] <- length(intersect(index_select, index_alter)) / length(index_alter)
          res_mat_list[[iter_method]][iter_simu, 2] <- length(setdiff(index_select, index_alter)) / length(index_select)
        }
      }
    }
    for (iter_method in seq_len(n_method)) {
      res_table_noconf[iter_method, 3] <-
        round(mean(res_mat_list[[iter_method]][, 1]), 4)
      res_table_noconf[iter_method, 4] <-
        round(sd(res_mat_list[[iter_method]][, 1]) / sqrt(n_simu), 4)
      res_table_noconf[iter_method, 5] <-
        round(mean(res_mat_list[[iter_method]][, 2]), 4)
      res_table_noconf[iter_method, 6] <-
        round(sd(res_mat_list[[iter_method]][, 2]) / sqrt(n_simu), 4)
    }

    #### save results
    res_mat[, (2 * (iter_ratio - 1) + 1):(2 * iter_ratio)] <- as.matrix(res_table_noconf[, c(3, 5)])
  }
  res_mat <- round(res_mat, 2)
  #### combine results
  cat(
    "LinDA & ", res_mat[1, 1], "&", res_mat[1, 2], "&", res_mat[1, 3], "&",
    res_mat[1, 4], "&", res_mat[1, 5], "&", res_mat[1, 6], "&", res_mat[1, 7], "&",
    res_mat[1, 8], "\\\\ \n"
  )
  cat(
    "LinDA97 & ", res_mat[2, 1], "&", res_mat[2, 2], "&", res_mat[2, 3], "&",
    res_mat[2, 4], "&", res_mat[2, 5], "&", res_mat[2, 6], "&", res_mat[2, 7], "&",
    res_mat[2, 8], "\\\\ \n"
  )
  cat(
    "LinDA90 & ", res_mat[3, 1], "&", res_mat[3, 2], "&", res_mat[3, 3], "&",
    res_mat[3, 4], "&", res_mat[3, 5], "&", res_mat[3, 6], "&", res_mat[3, 7], "&",
    res_mat[3, 8], "\\\\ \n"
  )
  cat(
    "Huber & ", res_mat[4, 1], "&", res_mat[4, 2], "&", res_mat[4, 3], "&",
    res_mat[4, 4], "&", res_mat[4, 5], "&", res_mat[4, 6], "&", res_mat[4, 7], "&",
    res_mat[4, 8], "\\\\ \n"
  )
  cat(
    "Bi-square & ", res_mat[5, 1], "&", res_mat[5, 2], "&", res_mat[5, 3], "&",
    res_mat[5, 4], "&", res_mat[5, 5], "&", res_mat[5, 6], "&", res_mat[5, 7], "&",
    res_mat[5, 8], "\\\\ \n"
  )
  cat(
    "QR & ", res_mat[6, 1], "&", res_mat[6, 2], "&", res_mat[6, 3], "&",
    res_mat[6, 4], "&", res_mat[6, 5], "&", res_mat[6, 6], "&", res_mat[6, 7], "&",
    res_mat[6, 8], "\\\\ \n"
  )
}
