set.seed(1)
library(ggplot2)
library(gridExtra)
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

####################### Run different methods #######################
formulas <- c("Disease", "Disease+Antibiotic", "Disease")
for (iter_dta in seq_len(3)) {
  ## load data set
  otu_tab <- dta[[2 * (iter_dta - 1) + 1]]
  meta <- dta[[2 * iter_dta]]
  dta_process <- preprocess_fun(otu_tab, meta)
  otu_tab <- dta_process$Y
  meta <- dta_process$Z
  formula <- formulas[iter_dta]
  out_res <- list()

  ## LinDA method
  linda_res <- LinDA::linda(otu_tab, meta, paste("~", formula))
  padj_linda <- linda_res$output[[1]]$padj
  out_res[["linda"]] <- padj_linda

  ## LinDA97
  ## remove samples for Y
  otu_tab_tmp <- robustDAA::winsor.fun(otu_tab, 0.97)
  N_tmp <- colSums(otu_tab_tmp)
  keep_sam <- which(N_tmp >= 1)
  otu_tab_use <- otu_tab[, keep_sam]
  ## remove samples for Z
  allvars <- colnames(meta)
  meta_use <- as.data.frame(meta[keep_sam, ])
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
  allvars <- colnames(meta)
  meta_use <- as.data.frame(meta[keep_sam, ])
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
    otu_tab = otu_tab, meta = meta, formula = paste("~", formula),
    reg_method = "psi.huber"
  )
  pvalue_huber <- huber_res$combine_pvalue
  padj_huber <- p.adjust(pvalue_huber, method = "BH")
  out_res[["huber"]] <- padj_huber

  ## Bisquare method
  bisquareres <- robustDAA::rlm_fun(
    otu_tab = otu_tab, meta = meta, formula = paste("~", formula),
    reg_method = "psi.bisquare"
  )
  pvalue_bisquare <- bisquareres$combine_pvalue
  padj_bisquare <- p.adjust(pvalue_bisquare, method = "BH")
  out_res[["bisquare"]] <- padj_bisquare

  ## QR method
  qrres <- robustDAA::qr_fun(otu_tab = otu_tab, meta = meta, paste("~", formula))
  pvalue_qr <- qrres$combine_pvalue
  padj_qr <- p.adjust(pvalue_qr, method = "BH")
  out_res[["qr"]] <- padj_qr

  ## save results
  saveRDS(out_res, paste0("real_data/results/data", iter_dta, "res.rds"))
}

####################### power plot #######################
thres_vec <- seq(0.01, 0.25, by = 0.01)
n_thres <- length(thres_vec)

## the first data set
dta_name1 <- names(dta)[1]
dta_name1 <- strsplit(dta_name1, "\\.")[[1]][1]
power_mat <- matrix(NA, nrow = 6, ncol = n_thres)
out_res <- readRDS(paste0("real_data/results/data", 1, "res.rds"))
padj_linda <- out_res[["linda"]]
padj_linda97 <- out_res[["linda97"]]
padj_linda90 <- out_res[["linda90"]]
padj_huber <- out_res[["huber"]]
padj_bisquare <- out_res[["bisquare"]]
padj_qr <- out_res[["qr"]]

for (iter_thres in seq_len(n_thres)) {
  thres_use <- thres_vec[iter_thres]
  power_mat[, iter_thres] <- c(
    sum(padj_linda <= thres_use),
    sum(padj_linda97 <= thres_use),
    sum(padj_linda90 <= thres_use),
    sum(padj_huber <= thres_use),
    sum(padj_bisquare <= thres_use),
    sum(padj_qr <= thres_use)
  )
}
df_plot1 <- data.frame(
  y = as.numeric(power_mat),
  method = rep(c("LinDA", "LinDA97", "LinDA90", "Huber", "Bi_square", "QR"), n_thres),
  x = rep(thres_vec, each = 6)
)
df_plot1$method <- factor(df_plot1$method, levels = c(
  "LinDA", "LinDA97", "LinDA90",
  "Huber", "Bi_square", "QR"
))
p1 <- ggplot(df_plot1, aes(x = x, y = y, color = method)) +
  geom_line(linewidth = 1, aes(linetype = method)) +
  geom_point(size = 1.5, aes(shape = method)) +
  xlab("Target FDR level") +
  ylab("Number of discoveries") +
  ggtitle(dta_name1) +
  theme_bw(base_size = 33) +
  theme(legend.position = "none")

## the second data set
dta_name2 <- names(dta)[3]
dta_name2 <- strsplit(dta_name2, "\\.")[[1]][1]
power_mat <- matrix(NA, nrow = 6, ncol = n_thres)
out_res <- readRDS(paste0("real_data/results/data", 2, "res.rds"))
padj_linda <- out_res[["linda"]]
padj_linda97 <- out_res[["linda97"]]
padj_linda90 <- out_res[["linda90"]]
padj_huber <- out_res[["huber"]]
padj_bisquare <- out_res[["bisquare"]]
padj_qr <- out_res[["qr"]]

for (iter_thres in seq_len(n_thres)) {
  thres_use <- thres_vec[iter_thres]
  power_mat[, iter_thres] <- c(
    sum(padj_linda <= thres_use),
    sum(padj_linda97 <= thres_use),
    sum(padj_linda90 <= thres_use),
    sum(padj_huber <= thres_use),
    sum(padj_bisquare <= thres_use),
    sum(padj_qr <= thres_use)
  )
}
df_plot2 <- data.frame(
  y = as.numeric(power_mat),
  method = rep(c("LinDA", "LinDA97", "LinDA90", "Huber", "Bi_square", "QR"), n_thres),
  x = rep(thres_vec, each = 6)
)
df_plot2$method <- factor(df_plot2$method, levels = c(
  "LinDA", "LinDA97", "LinDA90",
  "Huber", "Bi_square", "QR"
))
p2 <- ggplot(df_plot2, aes(x = x, y = y, color = method)) +
  geom_line(linewidth = 1, aes(linetype = method)) +
  geom_point(size = 1.5, aes(shape = method)) +
  xlab("Target FDR level") +
  ylab("Number of discoveries") +
  ggtitle(dta_name2) +
  theme_bw(base_size = 33) +
  theme(legend.position = "none")

## the third data set
dta_name3 <- names(dta)[5]
dta_name3 <- strsplit(dta_name3, "\\.")[[1]][1]
power_mat <- matrix(NA, nrow = 6, ncol = n_thres)
out_res <- readRDS(paste0("real_data/results/data", 3, "res.rds"))
padj_linda <- out_res[["linda"]]
padj_linda97 <- out_res[["linda97"]]
padj_linda90 <- out_res[["linda90"]]
padj_huber <- out_res[["huber"]]
padj_bisquare <- out_res[["bisquare"]]
padj_qr <- out_res[["qr"]]

for (iter_thres in seq_len(n_thres)) {
  thres_use <- thres_vec[iter_thres]
  power_mat[, iter_thres] <- c(
    sum(padj_linda <= thres_use),
    sum(padj_linda97 <= thres_use),
    sum(padj_linda90 <= thres_use),
    sum(padj_huber <= thres_use),
    sum(padj_bisquare <= thres_use),
    sum(padj_qr <= thres_use)
  )
}
df_plot3 <- data.frame(
  y = as.numeric(power_mat),
  method = rep(c("LinDA", "LinDA97", "LinDA90", "Huber", "Bi_square", "QR"), n_thres),
  x = rep(thres_vec, each = 6)
)
df_plot3$method <- factor(df_plot3$method, levels = c(
  "LinDA", "LinDA97", "LinDA90",
  "Huber", "Bi_square", "QR"
))
p3 <- ggplot(df_plot3, aes(x = x, y = y, color = method)) +
  geom_line(linewidth = 1, aes(linetype = method)) +
  geom_point(size = 1.5, aes(shape = method)) +
  xlab("Target FDR level") +
  ylab("Number of discoveries") +
  ggtitle(dta_name3) +
  theme_bw(base_size = 33) +
  theme(legend.position = "none")
p3_legend <- ggplot(df_plot3, aes(x = x, y = y, color = method)) +
  geom_line(linewidth = 1, aes(linetype = method)) +
  geom_point(size = 1.5, aes(shape = method)) +
  xlab("Target FDR level") +
  ylab("Number of discoveries") +
  ggtitle(dta_name3) +
  theme_bw(base_size = 33) +
  theme(legend.position = "bottom")
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}
shared_legend <- extract_legend(p3_legend)
## save plots
pdf("real_data/figures/power/res.pdf", width = 20, height = 10)
print(grid.arrange(arrangeGrob(p1, p2, p3, ncol = 3),
                   shared_legend,
                   nrow = 2, heights = c(10, 1)
))
dev.off()

####################### overlap plot #######################
library(UpSetR)
library(grid)
alpha <- 0.1
for (iter_dta in seq_len(3)) {
  dta_name <- names(dta)[2 * (iter_dta - 1) + 1]
  dta_name <- strsplit(dta_name, "\\.")[[1]][1]
  out_res <- readRDS(paste0("real_data/results/data", iter_dta, "res.rds"))
  padj_linda <- out_res[["linda"]]
  padj_linda97 <- out_res[["linda97"]]
  padj_linda90 <- out_res[["linda90"]]
  padj_huber <- out_res[["huber"]]
  padj_bisquare <- out_res[["bisquare"]]
  padj_qr <- out_res[["qr"]]
  res_list <- list()
  res_list[["LinDA"]] <- which(padj_linda < alpha)
  res_list[["LinDA97"]] <- which(padj_linda97 < alpha)
  res_list[["LinDA90"]] <- which(padj_linda90 < alpha)
  res_list[["Huber"]] <- which(padj_huber < alpha)
  res_list[["Bi_square"]] <- which(padj_bisquare < alpha)
  res_list[["QR"]] <- which(padj_qr < alpha)
  pdf(paste0("real_data/figures/overlap/realdata", iter_dta, ".pdf"),
    width = 10, height = 10
  )
  print(upset(fromList(res_list), nintersects = 12, mb.ratio = c(0.5, 0.5),
              sets.x.label = '# of Discoveries',order.by = c("freq", "degree"),
              decreasing = c(TRUE,FALSE), point.size = 3,
              text.scale = c(3, 3, 2, 3, 3, 3)))
  grid.text(dta_name, x = 0.65, y = 0.95, gp=gpar(fontsize = 22))
  dev.off()
}
