# =============================================================================
# 03_admixture_LEA.R
# sNMF admixture analysis for Cheirimedon femoratus
#
# Method: sNMF (sparse Non-negative Matrix Factorization)
#         implemented in the LEA R package
#         Reference: Frichot et al. (2014) Genetics 196:973-983
#
# Input:
#   final.recode.vcf  — filtered VCF
#   Amphi.txt         — population map
#
# Output:
#   output/cross_entropy_plot.pdf   — cross-entropy across K values
#   output/All_Qfiles/              — Q matrices for best runs
#   output/Figure_admixture_K*.pdf  — admixture bar plots per K
#
# Usage:
#   Rscript 03_admixture_LEA.R
#
# Requirements:
#   BiocManager::install("LEA")
#   install.packages("pophelper")
# =============================================================================

suppressPackageStartupMessages({
  library(LEA)
  library(pophelper)
  library(ggplot2)
})

dir.create("output/All_Qfiles", recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat(" sNMF admixture analysis — C. femoratus\n")
cat("============================================================\n")

# ── Parameters ────────────────────────────────────────────────────────────────
K_MAX       <- 10    # maximum K to test
REPS        <- 50    # repetitions per K value
CPU         <- 8     # parallel threads

# Color palette for admixture plots (up to K=10)
CLUSTER_COLS <- c(
  "#51692d", "#f1c039", "#f37d21", "#56ba32",
  "#a1d99b", "#9ecae1", "#fc9272", "#dd1c77",
  "#756bb1", "#636363"
)

# ── Step 1: Convert VCF to .geno format ──────────────────────────────────────

cat("\n[1/3] Converting VCF to .geno format...\n")

if (!file.exists("Qff.geno")) {
  vcf2geno(input.file  = "final.recode.vcf",
           output.file = "Qff.geno")
  cat("  Conversion complete: Qff.geno\n")
} else {
  cat("  Qff.geno already exists — skipping conversion.\n")
}

# ── Step 2: Run sNMF for K = 1 to K_MAX ──────────────────────────────────────
# entropy = TRUE: calculate cross-entropy for model selection
# project = "new": start fresh (use "continue" to add runs)

cat(sprintf("\n[2/3] Running sNMF (K=1 to %d, %d repetitions each)...\n",
            K_MAX, REPS))
cat("  This may take 30-60 minutes.\n\n")

project <- snmf(
  "Qff.geno",
  K           = 1:K_MAX,
  repetitions = REPS,
  entropy     = TRUE,
  CPU         = CPU,
  project     = "new"
)

# Cross-entropy plot
pdf("output/cross_entropy_plot.pdf", width = 6, height = 4)
plot(project, col = "maroon4", pch = 19, cex = 1.2,
     main = "Cross-entropy — C. femoratus sNMF",
     xlab = "K (number of ancestral populations)",
     ylab = "Cross-entropy")
dev.off()
cat("  Saved: output/cross_entropy_plot.pdf\n")

# Report best K per cross-entropy
cat("\n  Best run (lowest cross-entropy) per K:\n")
best_runs <- list()
for (k in 1:K_MAX) {
  best_k <- which.min(cross.entropy(project, K = k))
  best_runs[[k]] <- best_k
  cat(sprintf("    K=%d: run %d\n", k, best_k))
}

# ── Step 3: Export Q files for best runs and plot ─────────────────────────────

cat("\n[3/3] Exporting Q files and generating admixture plots...\n")

# Copy best Q files to All_Qfiles directory for pophelper
for (k in 2:K_MAX) {
  best_k  <- best_runs[[k]]
  q_file  <- Q(project, K = k, run = best_k)
  out_file <- sprintf("output/All_Qfiles/K%d_best_run%d.Q", k, best_k)
  write.table(q_file, out_file, col.names = FALSE, row.names = FALSE)
}

# Read all Q files with pophelper
sfiles <- list.files("output/All_Qfiles", full.names = TRUE, pattern = "\\.Q$")
sfiles <- sort(sfiles)
slist  <- readQ(files = sfiles)

# Generate admixture bar plots for K=2 to K=5
# Adjust K range and colours as needed for your data
k_plot_configs <- list(
  list(k=2, cols=CLUSTER_COLS[1:2]),
  list(k=3, cols=CLUSTER_COLS[1:3]),
  list(k=4, cols=CLUSTER_COLS[1:4]),
  list(k=5, cols=CLUSTER_COLS[1:5])
)

for (cfg in k_plot_configs) {
  k    <- cfg$k
  cols <- cfg$cols
  # Find the K=k file in slist
  idx  <- grep(paste0("K", k, "_"), names(slist))
  if (length(idx) == 0) next

  plotQ(
    qlist      = slist[idx],
    imgtype    = "pdf",
    height     = 1.5,
    clustercol = cols,
    dpi        = 300,
    exportpath = "output"
  )
  cat(sprintf("  Saved admixture plot: K=%d\n", k))
}

cat("\n============================================================\n")
cat(" sNMF analysis complete.\n")
cat(" Output files in: output/\n")
cat("   cross_entropy_plot.pdf   — model selection\n")
cat("   All_Qfiles/              — Q matrices for best runs\n")
cat("   Admixture bar plots per K\n")
cat("============================================================\n")
