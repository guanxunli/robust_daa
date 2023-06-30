library(GUniFrac)
library(parallel)
n_simu <- 100
## load reference data
data(stool.otu.tab)
comm <- stool.otu.tab

## define parameters
n_sam <- 50
n_taxa <- 500
set_df <- data.frame(n_sam = rep(n_sam, length(n_taxa)), n_taxa = rep(n_taxa, each = length(n_sam)))
signa_den <- c(0.05, 0.2)
set_df <- cbind(apply(set_df, 2, rep, length(signa_den)), rep(signa_den, each = nrow(set_df)))
signa_streng <- c(1, 1.25, 1.5, 1.75, 2)
set_df <- cbind(apply(set_df, 2, rep, length(signa_streng)), rep(signa_streng, each = nrow(set_df)))
diff_mode <- c("abundant", "mix")
set_df <- cbind(apply(set_df, 2, rep, length(diff_mode)), rep(diff_mode, each = nrow(set_df)))
colnames(set_df) <- c("n_sam", "n_taxa", "signa_den", "signa_streng", "mode")
set_df <- as.data.frame(set_df)
nset <- nrow(set_df)

#### without confounder
for (iter in seq_len(nset)) {
  print(iter)
  n_sam <- as.numeric(set_df$n_sam[iter])
  n_taxa <- as.numeric(set_df$n_taxa[iter])
  signa_den <- as.numeric(set_df$signa_den[iter])
  signa_streng <- as.numeric(set_df$signa_streng[iter])
  diff_mode <- set_df$mode[iter]
  set.seed(1)
  dta_list <- mclapply(seq_len(n_simu), function(x) {
    sim.obj <- SimulateMSeq(
      ref.otu.tab = comm, nSam = n_sam, nOTU = n_taxa,
      diff.otu.pct = signa_den, confounder.type = "none",
      covariate.eff.mean = signa_streng, diff.otu.mode = diff_mode
    )
  }, mc.cores = 50)
  saveRDS(dta_list, paste0(
    "stool_data/outlier0/datasets/noconf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, "mode",  diff_mode, ".rds"
  ))
}