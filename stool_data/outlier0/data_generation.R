set.seed(1)
library(GUniFrac)
library(parallel)
n_simu <- 100
## load reference data
data(stool.otu.tab)
comm <- stool.otu.tab

## define parameters
n_sam <- c(50, 100, 200)
n_taxa <- c(50, 100, 500)
set_df <- data.frame(n_sam = rep(n_sam, length(n_taxa)), n_taxa = rep(n_taxa, each = length(n_sam)))
signa_den <- c(0.05, 0.2)
set_df <- cbind(apply(set_df, 2, rep, length(signa_den)), rep(signa_den, each = nrow(set_df)))
signa_streng <- 1.25
set_df <- cbind(apply(set_df, 2, rep, length(signa_streng)), rep(signa_streng, each = nrow(set_df)))
colnames(set_df) <- c("n_sam", "n_taxa", "signa_den", "signa_streng")
set_df <- as.data.frame(set_df)
nset <- nrow(set_df)

#### without confounder
for (iter in seq_len(nset)) {
  n_sam <- as.numeric(set_df$n_sam[iter])
  n_taxa <- as.numeric(set_df$n_taxa[iter])
  signa_den <- as.numeric(set_df$signa_den[iter])
  signa_streng <- as.numeric(set_df$signa_streng[iter])
  dta_list <- mclapply(seq_len(n_simu), function(x) {
    sim.obj <- SimulateMSeq(
      ref.otu.tab = comm, nSam = n_sam, nOTU = n_taxa,
      diff.otu.pct = signa_den, confounder.type = "none",
      confounder.eff.mean = signa_streng
    )
  }, mc.cores = 50)
  saveRDS(dta_list, paste0(
    "stool_data/outlier0/datasets/noconf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))
}

#### with confounder
for (iter in seq_len(nset)) {
  n_sam <- as.numeric(set_df$n_sam[iter])
  n_taxa <- as.numeric(set_df$n_taxa[iter])
  signa_den <- as.numeric(set_df$signa_den[iter])
  signa_streng <- as.numeric(set_df$signa_streng[iter])
  dta_list <- mclapply(seq_len(n_simu), function(x) {
    sim.obj <- SimulateMSeq(
      ref.otu.tab = comm, nSam = n_sam, nOTU = n_taxa,
      diff.otu.pct = signa_den, confounder.type = "both",
      confounder.eff.mean = signa_streng,
      conf.diff.otu.pct = 1, conf.nondiff.otu.pct = 0.1
    )
  }, mc.cores = 50)
  saveRDS(dta_list, paste0(
    "stool_data/outlier0/datasets/conf_nsam", n_sam, "ntaxa", n_taxa,
    "signal", signa_den, "streng", signa_streng, ".rds"
  ))
}
