set.seed(1)
################## define parameters ##################
library(ggplot2)
library(gridExtra)
method_vec <- c(
  "linda", "linda97", "linda90", "linda_winsor",
  "huber", "bisquare"
)
n_method <- length(method_vec)
outlier <- "outlier2"
## load parameters
para0 <- readRDS(paste0("loglinear_data/", outlier, "/datasets/log.normal.para.rds"))
beta0 <- para0$beta0
sigma2 <- para0$sigma2
# parameter use
sample_size_vec <- c(50, 50, 200, 200)
m <- 500
n_sim <- 100
# define settings
sig_density_vec <- c(0.05, 0.2, 0.05, 0.2)
sig_strength_vec <- seq(1.05, 2, length.out = 6)
n_signa <- length(sig_strength_vec)
conf_vec <- c("nocon_n", "con_n")

################## without confounder ##################
for (iter_conf in conf_vec) {
  for (iter_plot in seq_len(4)) {
    n <- sample_size_vec[iter_plot]
    gamma <- sig_density_vec[iter_plot]
    out_res_list <- list()
    for (iter_sig in seq_len(n_signa)) {
      out_res_list[[iter_sig]] <- list()
      mu_use <- sig_strength_vec[iter_sig]
      ## load results
      dta_list <- readRDS(paste0(
        "loglinear_data/", outlier, "/datasets/", iter_conf, n,
        "gamma", gamma, "mu", mu_use, ".rds"
      ))
      res_list <- list()
      res_mat_list <- list()
      for (iter_method in method_vec) {
        res_list[[iter_method]] <- readRDS(paste0(
          "loglinear_data/", outlier, "/results/", iter_method, "_",
          iter_conf, n, "gamma", gamma, "mu", mu_use, ".rds"
        ))
        res_mat_list[[iter_method]] <- matrix(NA, nrow = n_sim, ncol = 2)
        colnames(res_mat_list[[iter_method]]) <- c("power", "fdr")
      }
      ## calculate results
      for (iter_simu in seq_len(n_sim)) {
        dta <- dta_list[[iter_simu]]
        index_alter <- dta$index_alter
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
      out_res_list[[iter_sig]] <- res_mat_list
    }
    ## rewrite results
    power_list <- list()
    power_sd_list <- list()
    fdr_list <- list()
    fdr_sd_list <- list()
    for (iter_method in method_vec) {
      power_list[[iter_method]] <- numeric(n_signa)
      fdr_list[[iter_method]] <- numeric(n_signa)
      power_sd_list[[iter_method]] <- numeric(n_signa)
      fdr_sd_list[[iter_method]] <- numeric(n_signa)
      for (iter_sig in seq_len(n_signa)) {
        power_list[[iter_method]][iter_sig] <- mean(out_res_list[[iter_sig]][[iter_method]][, 1])
        power_sd_list[[iter_method]][iter_sig] <- sd(out_res_list[[iter_sig]][[iter_method]][, 1])
        fdr_list[[iter_method]][iter_sig] <- mean(out_res_list[[iter_sig]][[iter_method]][, 2])
        fdr_sd_list[[iter_method]][iter_sig] <- mean(out_res_list[[iter_sig]][[iter_method]][, 2])
      }
    }
    ## data frame to plots
    df_plot <- data.frame(
      power = unlist(power_list), fdr = unlist(fdr_list), power_sd = unlist(power_sd_list),
      fdr_sd = unlist(fdr_sd_list), sig_streng = rep(sig_strength_vec, length(method_vec)),
      method = rep(method_vec, each = n_signa)
    )
    df_plot$method[which(df_plot$method == "linda")] <- "LinDA"
    df_plot$method[which(df_plot$method == "linda97")] <- "LinDA97"
    df_plot$method[which(df_plot$method == "linda90")] <- "LinDA90"
    df_plot$method[which(df_plot$method == "linda_winsor")] <- "LinDA_winsor"
    df_plot$method[which(df_plot$method == "huber")] <- "Huber"
    df_plot$method[which(df_plot$method == "bisquare")] <- "Bi_square"
    # df_plot$method[which(df_plot$method == "qr")] <- "QR"
    df_plot$method <- factor(df_plot$method, levels = c(
      "LinDA", "LinDA97", "LinDA90", "LinDA_winsor",
      "Huber", "Bi_square"
    ))
    p1 <- ggplot(data = df_plot, aes(x = sig_streng, y = power, color = method)) +
      geom_line(aes(linetype = method), linewidth = 1.5) +
      geom_point() +
      xlab("Signal strength") +
      ylab("Power") +
      xlim(c(1, 2)) +
      theme_bw(base_size = 22) +
      theme(legend.position = "none")
    p2 <- ggplot(data = df_plot, aes(x = sig_streng, y = fdr, color = method)) +
      geom_line(aes(linetype = method), linewidth = 1.5) +
      geom_point() +
      xlab("Signal strength") +
      xlim(c(1, 2)) +
      ylab("FDR") +
      geom_hline(yintercept = 0.05, linetype = "dashed") +
      theme_bw(base_size = 22) +
      theme(legend.position = "none")
    p2_legends <- ggplot(data = df_plot, aes(x = sig_streng, y = fdr, color = method)) +
      geom_line(aes(linetype = method), linewidth = 1.5) +
      geom_point() +
      xlab("Signal strength") +
      xlim(c(1, 2)) +
      ylab("FDR") +
      geom_hline(yintercept = 0.05, linetype = "dashed") +
      theme_bw(base_size = 22) +
      theme(legend.position = "bottom")
    extract_legend <- function(my_ggp) {
      step1 <- ggplot_gtable(ggplot_build(my_ggp))
      step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
      step3 <- step1$grobs[[step2]]
      return(step3)
    }
    shared_legend <- extract_legend(p2_legends)
    ## save plots
    pdf(paste0(
      "loglinear_data/", outlier, "/figures/", iter_conf, n,
      "gamma", gamma, ".pdf"
    ), width = 10, height = 10)
    print(grid.arrange(arrangeGrob(p1, p2, nrow = 2),
      shared_legend,
      nrow = 2, heights = c(10, 1)
    ))
    dev.off()
  }
}
