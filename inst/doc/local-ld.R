## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 4)
## This vignette intentionally uses tiny marker counts to keep its examples
## fast and to illustrate finite-locus behavior; the warning is explained in
## the text below instead of repeated in every chunk's output.
options(rBahadur.warn_small_m = FALSE)

## -----------------------------------------------------------------------------
library(rBahadur)

panel <- kg_reference()
str(panel[c("pos", "cM", "chrom", "build")], max.level = 1)
dim(panel$haplotypes)

## -----------------------------------------------------------------------------
H <- matrix(as.integer(panel$haplotypes), nrow = nrow(panel$haplotypes))

ld_decay <- function(G, pos, breaks = c(0, 5, 10, 25, 50, 100, 250, 1000)) {
  set.seed(1)
  ii <- sort(sample(ncol(G), 300))
  R <- suppressWarnings(stats::cor(G[, ii]))
  d <- abs(outer(pos[ii], pos[ii], "-")) / 1000
  ut <- upper.tri(R)
  tapply(R[ut]^2, cut(d[ut], breaks), mean, na.rm = TRUE)
}

round(ld_decay(H, panel$pos), 4)

## -----------------------------------------------------------------------------
set.seed(2026)
sim <- am_mosaic(h2_0 = 0.5, r = 0.5, n = 1500, panel = panel, m = 10)

dim(sim$X)
length(sim$causal_idx)

## -----------------------------------------------------------------------------
g_from_X <- as.vector(sim$X[, sim$causal_idx, drop = FALSE] %*% sim$beta_raw)
all.equal(g_from_X, as.vector(sim$g))

## -----------------------------------------------------------------------------
c(empirical = var(as.vector(sim$g)), infinitesimal_limit = vg_eq(0.5, 0.5, 0.5))

## -----------------------------------------------------------------------------
set.seed(11)
big <- am_mosaic(h2_0 = 0.5, r = 0.5, n = 1500, panel = panel, m = 200)
c(empirical = var(as.vector(big$g)), infinitesimal_limit = vg_eq(0.5, 0.5, 0.5))
c(empirical = var(as.vector(big$g)) / var(as.vector(big$y)),
  infinitesimal_limit = h2_eq(0.5, 0.5))

## -----------------------------------------------------------------------------
round(ld_decay(sim$X, panel$pos), 4)

## -----------------------------------------------------------------------------
set.seed(2026)
flat <- am_simulate(h2_0 = 0.5, r = 0.5, m = 500, n = 1500)
round(ld_decay(flat$X, panel$pos[seq_len(500)]), 4)

## -----------------------------------------------------------------------------
reach <- sapply(c(5, 20, 80), function(m) {
  set.seed(3)
  s <- am_mosaic(0.5, 0.3, n = 400, panel = panel, m = m)
  d <- ld_decay(s$X, panel$pos)
  d[["(10,25]"]]
})
data.frame(m = c(5, 20, 80),
           mean_block_kb = round(diff(range(panel$pos)) / 1000 / c(5, 20, 80)),
           r2_at_10_25kb = round(reach, 4))

## ----eval = FALSE-------------------------------------------------------------
#  ## build a panel directly from files you already have
#  panel <- vcf_to_panel(
#    vcf = "ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz",
#    map = "plink.chr22.GRCh38.map",
#    min_maf = 0.01
#  )
#  
#  ## or fetch a region of 1000 Genomes and its genetic map in one call
#  panel <- download_1kg_panel(chrom = "22", start = 20e6, end = 30e6)

## ----eval = FALSE-------------------------------------------------------------
#  out <- am_mosaic(h2_0 = 0.5, r = 0.5, n = 50000, panel = panel, m = 2000,
#                   path = "chr22_am", format = "bed")
#  
#  ## readable by plink, GCTA, or any other standard tool
#  ## and by read_genotypes() back in R
#  ## the .bim preserves chromosome, position, map, ID, and allele metadata

