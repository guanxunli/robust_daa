set.seed(1)
outlier <- "outlier0"
################## define parameters ##################
library(ggplot2)
library(gridExtra)
method_vec <- c(
  "linda", "linda97", "linda90", "linda_winsor", "huber", "bisquare", "qr"
)
n_method <- length(method_vec)
n_sam <- c(50, 100, 200)
n_taxa <- c(50, 100, 500)
set_df <- data.frame(n_sam = rep(n_sam, length(n_taxa)), n_taxa = rep(n_taxa, each = length(n_sam)))
signa_den <- c(0.05, 0.2)
set_df <- cbind(apply(set_df, 2, rep, length(signa_den)), rep(signa_den, each = nrow(set_df)))
nset <- nrow(set_df)
set_df <- as.data.frame(set_df)
colnames(set_df) <- c("n_sam", "n_taxa", "signa_den")
signa_streng <- 1.25
n_simu <- 1e2
################## without confounder ##################
for (iter_set in seq_len(nset)) {
  print(iter_set)
  n_sam <- as.numeric(set_df$n_sam[iter_set])
  n_taxa <- as.numeric(set_df$n_taxa[iter_set])
  signa_den <- as.numeric(set_df$signa_den[iter_set])
  #### without confounders
  dta_list <- readRDS(paste0(
    "stool_data/", outlier, "/datasets/noconf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))
  res_table_noconf <- data.frame(
    "method" = c("LinDA", "LinDA97", "LinDA90", "LinDA_winsor", "Huber", "Bi-square", "QR"),
    "confounder" = rep("without", n_method),
    "Power" = numeric(n_method), "Power_sd" = numeric(n_method),
    "FDR" = numeric(n_method), "FDR_sd" = numeric(n_method)
  )
  ## load results
  res_list <- list()
  res_mat_list <- list()
  for (iter_method in method_vec) {
    res_list[[iter_method]] <- readRDS(paste0(
      "stool_data/", outlier, "/results/", iter_method, "_noconf_nsam", n_sam, "ntaxa", n_taxa,
      "signal", signa_den, "streng", signa_streng, ".rds"
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

  #### with confounders
  dta_list <- readRDS(paste0(
    "stool_data/", outlier, "/datasets/conf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))
  res_table_conf <- data.frame(
    "method" = c("LinDA", "LinDA97", "LinDA90", "LinDA_winsor", "Huber", "Bi-square", "QR"),
    "confounder" = rep("with", n_method),
    "Power" = numeric(n_method), "Power_sd" = numeric(n_method),
    "FDR" = numeric(n_method), "FDR_sd" = numeric(n_method)
  )
  ## load results
  res_list <- list()
  res_mat_list <- list()
  for (iter_method in method_vec) {
    res_list[[iter_method]] <- readRDS(paste0(
      "stool_data/", outlier, "/results/", iter_method, "_conf_nsam", n_sam, "ntaxa", n_taxa,
      "signal", signa_den, "streng", signa_streng, ".rds"
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
    res_table_conf[iter_method, 3] <-
      round(mean(res_mat_list[[iter_method]][, 1]), 4)
    res_table_conf[iter_method, 4] <-
      round(sd(res_mat_list[[iter_method]][, 1]) / sqrt(n_simu), 4)
    res_table_conf[iter_method, 5] <-
      round(mean(res_mat_list[[iter_method]][, 2]), 4)
    res_table_conf[iter_method, 6] <-
      round(sd(res_mat_list[[iter_method]][, 2]) / sqrt(n_simu), 4)
  }

  #### combine results
  res_table <- rbind(res_table_noconf, res_table_conf)
  res_table <- insight::export_table(res_table, format = "md")
  write.table(res_table, file = paste0(
    "stool_data/", outlier,
    "/tables/res_", n_sam, "_", n_taxa, "_", signa_den, ".txt"
  ))
}
