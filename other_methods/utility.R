library(ALDEx2)
library(ANCOMBC)
library(phyloseq)
library(DESeq2)
library(GMPR)
library(Maaslin2)

aldex2.fun <- function(otu.tab, meta, formula, alpha = 0.05) {
  tt <- FALSE
  if (length(all.vars(as.formula(paste("~", formula)))) == 1) {
    if (!is.numeric(meta[, formula])) {
      tt <- TRUE
    }
  }
  if (tt) {
    res <- aldex(otu.tab, meta[, formula], test = "t", effect = FALSE, denom = "all")
    pval <- res$wi.ep
  } else {
    design <- model.matrix(as.formula(paste("~", formula)), meta)
    x <- aldex.clr(otu.tab, design, denom = "all")
    res <- aldex.glm(x, design)
    pval <- res[, 8]
  }
  qval <- p.adjust(pval, method = "BH")
  rej <- which(qval <= alpha)
  if (length(pval) != nrow(otu.tab)) {
    rej <- as.numeric(gsub("taxon", "", rownames(res)))[rej]
  }
  return(rej)
}

## If use **zero_cut = 1**, the taxa with 100% zeros will still be removed.
## So use **zero_cut = 1.1** to ensure to keep all the taxa.
## Use **lib_cut = 1** to keep all the samples.
ancombc.fun <- function(otu.tab, meta, formula, alpha = 0.05) {
  OTU <- otu_table(otu.tab, taxa_are_rows = TRUE)
  META <- sample_data(meta)
  PHYSEQ <- phyloseq(OTU, META)
  
  res <- ancombc(
    phyloseq = PHYSEQ, formula = formula, p_adj_method = "BH",
    zero_cut = 1.1, lib_cut = 1, conserve = TRUE
  )
  
  pval <- res$res$p_val[, 1]
  qval <- p.adjust(pval, method = "BH")
  rej <- which(qval <= alpha)
  return(rej)
}

maaslin2.fun <- function(otu.tab, meta, fix.eff, rand.eff = NULL, alpha = 0.05) {
  res <- Maaslin2(
    input_data = otu.tab,
    input_metadata = meta,
    output = "masslin2",
    min_prevalence = 0,
    fixed_effects = fix.eff,
    random_effects = rand.eff,
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )
  ind <- which(res$results$metadata == fix.eff[1])
  ord <- match(rownames(otu.tab), res$results$feature[ind])
  pval <- res$results$pval[ind][ord]
  qval <- p.adjust(pval, method = "BH")
  rej <- which(qval <= alpha)
  return(rej)
}