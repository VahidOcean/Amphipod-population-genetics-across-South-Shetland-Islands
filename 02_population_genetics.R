# =============================================================================
# 02_population_genetics.R
# Population genetic analyses for Cheirimedon femoratus
# South Shetland Islands — four populations (DICF, KGCF, LICF, SICF)
#
# Analyses:
#   1.  PCA (principal components analysis)
#   2.  FST (pairwise differentiation)
#   3.  AMOVA (analysis of molecular variance)
#   4.  Admixture / sNMF clustering (LEA)
#   5.  Observed and expected heterozygosity (He, Ho)
#   6.  UPGMA phylogenetic tree
#   7.  Minimum spanning network
#   8.  DAPC (discriminant analysis of principal components)
#   9.  PCoA (principal coordinates analysis)
#
# Input files:
#   final.recode.vcf  — filtered VCF (from 00_stacks_vcftools_pipeline.sh)
#   Amphi.txt         — population map (sample_id  population)
#
# Output:
#   Figures and tables saved to ./output/
#
# Requirements:
#   install.packages(c("vcfR","adegenet","hierfstat","poppr","StAMPP",
#                      "ggplot2","reshape2","vegan","ape","RColorBrewer",
#                      "ggrepel","tidyr","LEA","pophelper"))
#
# Usage:
#   Rscript 02_population_genetics.R
#   Or run interactively in RStudio
# =============================================================================

# ── 0. Setup ──────────────────────────────────────────────────────────────────

# Set working directory if running interactively
# setwd("/path/to/your/data")

# Create output directory
dir.create("output", showWarnings = FALSE)

# Load required packages
suppressPackageStartupMessages({
  library(vcfR)
  library(adegenet)
  library(hierfstat)
  library(poppr)
  library(StAMPP)
  library(ggplot2)
  library(reshape2)
  library(ape)
  library(RColorBrewer)
  library(ggrepel)
  library(tidyr)
})

# Population color palette — consistent across all figures
POP_COLORS <- c(
  DICF = "#f1c039",
  KGCF = "#f37d21",
  LICF = "#51692d",
  SICF = "#56ba32"
)

POPULATIONS <- c("DICF", "KGCF", "LICF", "SICF")

cat("============================================================\n")
cat(" Population genetics analysis — C. femoratus\n")
cat("============================================================\n")


# ── 1. Load data ──────────────────────────────────────────────────────────────

cat("\n[1/9] Loading VCF and population map...\n")

vcf      <- read.vcfR("final.recode.vcf", verbose = FALSE)
pop_data <- read.table("Amphi.txt", header = FALSE,
                       col.names = c("SampleID", "Population"),
                       stringsAsFactors = FALSE)

# Convert VCF to genlight (for PCA, FST, tree, MSN)
gl <- vcfR2genlight(vcf)
ploidy(gl) <- 2
pop(gl)    <- pop_data$Population[match(indNames(gl), pop_data$SampleID)]

# Convert VCF to genind (for AMOVA, He/Ho, DAPC)
genind_obj <- vcfR2genind(vcf)
pop(genind_obj) <- pop_data$Population[match(indNames(genind_obj),
                                              pop_data$SampleID)]

cat(sprintf("  Individuals : %d\n", nInd(gl)))
cat(sprintf("  SNPs        : %d\n", nLoc(gl)))
cat(sprintf("  Populations : %s\n", paste(levels(pop(gl)), collapse=", ")))


# ── 2. Principal Components Analysis (PCA) ────────────────────────────────────
# Method: glPca from adegenet
# Reference: Jombart & Ahmed (2011) Bioinformatics 27:3070-3071

cat("\n[2/9] Running PCA...\n")

pca <- glPca(gl, nf = 10)

# Variance explained per PC
pct_var <- round(100 * pca$eig / sum(pca$eig), 2)
cat(sprintf("  PC1: %.2f%%  PC2: %.2f%%\n", pct_var[1], pct_var[2]))

# Save eigenvalues
write.csv(
  data.frame(PC = seq_along(pct_var), Eigenvalue = pca$eig,
             Percent_variance = pct_var),
  "output/PCA_eigenvalues.csv", row.names = FALSE
)

# PCA scores dataframe
pca_scores <- as.data.frame(pca$scores)
pca_scores$Population <- pop(gl)

# Plot
p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 2.5, alpha = 0.85) +
  stat_ellipse(level = 0.95, linewidth = 0.7, linetype = 2) +
  scale_color_manual(values = POP_COLORS) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
  labs(
    x = paste0("PC1 (", pct_var[1], "%)"),
    y = paste0("PC2 (", pct_var[2], "%)"),
    color = "Population"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

ggsave("output/Figure_PCA.pdf", p_pca, width = 6, height = 4.5)
ggsave("output/Figure_PCA.png", p_pca, width = 6, height = 4.5, dpi = 300)
cat("  Saved: output/Figure_PCA.pdf\n")


# ── 3. Pairwise FST ───────────────────────────────────────────────────────────
# Method: StAMPP (Pembleton et al. 2013)
# 95% CI via bootstrapping (100 replicates)
# Reference: Weir & Cockerham (1984) FST estimator

cat("\n[3/9] Calculating pairwise FST...\n")

fst_result <- stamppFst(gl, nboots = 100, percent = 95, nclusters = 1)
fst_matrix  <- fst_result$Fsts
fst_pvalues <- fst_result$Pvalues

write.table(fst_matrix,  "output/FST_matrix.txt",  sep = "\t", quote = FALSE)
write.table(fst_pvalues, "output/FST_pvalues.txt", sep = "\t", quote = FALSE)

cat("  Pairwise FST matrix:\n")
print(round(fst_matrix, 4))
cat("  P-values:\n")
print(fst_pvalues)

# FST heatmap
fst_melt <- melt(as.matrix(fst_matrix), na.rm = TRUE)

p_fst <- ggplot(fst_melt, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(value, 3)), size = 3.5,
            color = "white", fontface = "bold") +
  scale_fill_gradient2(
    low     = "#ffd60a",
    mid     = "#4e9de6",
    high    = "#001d3d",
    midpoint = 0.013,
    name    = expression(italic(F)[ST])
  ) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

ggsave("output/Figure_FST_heatmap.pdf", p_fst, width = 5, height = 4.5)
ggsave("output/Figure_FST_heatmap.png", p_fst, width = 5, height = 4.5, dpi = 300)
cat("  Saved: output/Figure_FST_heatmap.pdf\n")


# ── 4. AMOVA ─────────────────────────────────────────────────────────────────
# Analysis of Molecular Variance
# Method: poppr.amova (Excoffier et al. 1992)
# Significance: 1000 random permutations

cat("\n[4/9] Running AMOVA...\n")

amova_result <- poppr.amova(genind_obj, ~Population, within = FALSE)
cat("  Variance components:\n")
print(amova_result$componentsofcovariance)

set.seed(123)
amova_test <- randtest(amova_result, nrepet = 1000)
cat(sprintf("  AMOVA p-value: %.4f\n", amova_test$pvalue))

# Save results
sink("output/AMOVA_results.txt")
cat("AMOVA results — C. femoratus\n")
cat("==============================\n")
print(amova_result$componentsofcovariance)
cat("\nSignificance test (1000 permutations):\n")
print(amova_test)
sink()
cat("  Saved: output/AMOVA_results.txt\n")

# AMOVA variance components bar chart
var_pct <- amova_result$componentsofcovariance$`%`
var_df  <- data.frame(
  Component  = c("Among populations", "Within populations"),
  Percentage = var_pct[1:2]
)

p_amova <- ggplot(var_df, aes(x = Component, y = Percentage, fill = Component)) +
  geom_bar(stat = "identity", width = 0.5, color = "black", linewidth = 0.4) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("#2166ac", "#d1e5f0")) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(y = "Variance explained (%)", x = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

ggsave("output/Figure_AMOVA.pdf", p_amova, width = 4.5, height = 4)
ggsave("output/Figure_AMOVA.png", p_amova, width = 4.5, height = 4, dpi = 300)
cat("  Saved: output/Figure_AMOVA.pdf\n")


# ── 5. Heterozygosity (He and Ho) ─────────────────────────────────────────────
# Method: basic.stats from hierfstat
# He = expected heterozygosity; Ho = observed heterozygosity
# FIS = (He - Ho) / He

cat("\n[5/9] Calculating heterozygosity (He, Ho, FIS)...\n")

hierf_data  <- genind2hierfstat(genind_obj)
basic       <- basic.stats(hierf_data, diploid = TRUE)

Ho_means <- colMeans(basic$Ho, na.rm = TRUE)
He_means <- colMeans(basic$Hs, na.rm = TRUE)
Fis_vals <- (He_means - Ho_means) / He_means

het_df <- data.frame(
  Population = names(Ho_means),
  Ho         = round(Ho_means, 4),
  He         = round(He_means, 4),
  FIS        = round(Fis_vals, 4)
)
het_df <- het_df[match(POPULATIONS, het_df$Population), ]

cat("  Heterozygosity summary:\n")
print(het_df)
write.csv(het_df, "output/Heterozygosity_summary.csv", row.names = FALSE)

# He vs Ho bar chart
het_melt <- melt(het_df[, c("Population","Ho","He")], id.vars = "Population")
het_melt$Population <- factor(het_melt$Population, levels = POPULATIONS)

p_het <- ggplot(het_melt, aes(x = Population, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.65, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c(Ho = "#4E79A7", He = "#F28E2B"),
                    labels = c(Ho = "Observed (Ho)", He = "Expected (He)"),
                    name = NULL) +
  labs(y = "Heterozygosity", x = "Population") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top")

ggsave("output/Figure_Heterozygosity.pdf", p_het, width = 5.5, height = 4)
ggsave("output/Figure_Heterozygosity.png", p_het, width = 5.5, height = 4, dpi = 300)
cat("  Saved: output/Figure_Heterozygosity.pdf\n")


# ── 6. UPGMA phylogenetic tree ────────────────────────────────────────────────
# Method: aboot (poppr) using bitwise distance, 100 bootstrap replicates
# Tree type: UPGMA

cat("\n[6/9] Building UPGMA bootstrap tree...\n")

set.seed(42)
tree <- aboot(gl, tree = "upgma", distance = bitwise.dist,
              sample = 100, showtree = FALSE, cutoff = 50, quiet = TRUE)

# Tip labels: sample - population
tree$tip.label <- paste(indNames(gl), pop(gl), sep = " — ")

# Save Newick tree
write.tree(tree, "output/UPGMA_tree.nwk")

# Plot
pdf("output/Figure_UPGMA_tree.pdf", width = 8, height = 10)
n_pop  <- nlevels(pop(gl))
cols   <- brewer.pal(max(3, n_pop), "Dark2")[seq_len(n_pop)]
names(cols) <- levels(pop(gl))
plot.phylo(tree, cex = 0.6, font = 2, adj = 0,
           tip.color = cols[pop(gl)])
nodelabels(tree$node.label, adj = c(1.3, -0.5),
           frame = "n", cex = 0.6, font = 3, xpd = TRUE)
legend("topleft", legend = levels(pop(gl)),
       fill = cols, border = FALSE, bty = "n", cex = 0.9)
axis(side = 1)
title(xlab = "Genetic distance (proportion of different loci)")
dev.off()
cat("  Saved: output/Figure_UPGMA_tree.pdf\n")


# ── 7. Minimum spanning network ───────────────────────────────────────────────

cat("\n[7/9] Building minimum spanning network...\n")

rubi_dist <- bitwise.dist(gl)
rubi_msn  <- poppr.msn(gl, rubi_dist, showplot = FALSE, include.ties = TRUE)

# Set node sizes
node_size <- rep(2, nInd(gl))
names(node_size) <- indNames(gl)
vertex.attributes(rubi_msn$graph)$size <- node_size

pdf("output/Figure_MSN.pdf", width = 7, height = 7)
set.seed(42)
plot_poppr_msn(gl, rubi_msn,
               palette = cols,
               gadj    = 70)
dev.off()
cat("  Saved: output/Figure_MSN.pdf\n")


# ── 8. DAPC ───────────────────────────────────────────────────────────────────
# Discriminant Analysis of Principal Components
# Reference: Jombart et al. (2010) BMC Genetics

cat("\n[8/9] Running DAPC...\n")

set.seed(42)
dapc_result <- dapc(gl, n.pca = 10, n.da = 3)

# Posterior membership probabilities
dapc_df <- as.data.frame(dapc_result$posterior)
dapc_df$Population  <- pop(gl)
dapc_df$Sample      <- rownames(dapc_df)
dapc_long <- pivot_longer(dapc_df, -c(Population, Sample),
                           names_to = "Assigned_Pop",
                           values_to = "Posterior_prob")

p_dapc <- ggplot(dapc_long,
                 aes(x = Sample, y = Posterior_prob, fill = Assigned_Pop)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = POP_COLORS, name = "Assigned population") +
  facet_grid(~Population, scales = "free_x", space = "free_x") +
  labs(y = "Posterior membership probability", x = NULL) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1, size = 5),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    legend.position  = "top"
  )

ggsave("output/Figure_DAPC.pdf", p_dapc, width = 9, height = 3.5)
ggsave("output/Figure_DAPC.png", p_dapc, width = 9, height = 3.5, dpi = 300)
cat("  Saved: output/Figure_DAPC.pdf\n")


# ── 9. PCoA ───────────────────────────────────────────────────────────────────
# Principal Coordinates Analysis on individual-level genetic distances

cat("\n[9/9] Running PCoA...\n")

gen_dist    <- dist(tab(genind_obj, NA.method = "mean"))
pcoa_result <- cmdscale(gen_dist, k = 2, eig = TRUE)
pco_var     <- round(pcoa_result$eig / sum(abs(pcoa_result$eig)) * 100, 2)[1:2]

pcoa_df <- data.frame(
  Sample     = rownames(pcoa_result$points),
  PCo1       = pcoa_result$points[, 1],
  PCo2       = pcoa_result$points[, 2],
  Population = pop(genind_obj)
)

p_pcoa <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2, color = Population)) +
  geom_point(size = 2.5, alpha = 0.85) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.6) +
  scale_color_manual(values = POP_COLORS) +
  labs(
    x     = paste0("PCo1 (", pco_var[1], "%)"),
    y     = paste0("PCo2 (", pco_var[2], "%)"),
    color = "Population"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    panel.border    = element_rect(fill = NA, color = "black", linewidth = 0.5)
  )

ggsave("output/Figure_PCoA.pdf", p_pcoa, width = 6, height = 4.5)
ggsave("output/Figure_PCoA.png", p_pcoa, width = 6, height = 4.5, dpi = 300)
cat("  Saved: output/Figure_PCoA.pdf\n")


# ── Summary ───────────────────────────────────────────────────────────────────

cat("\n============================================================\n")
cat(" All analyses complete. Output files in: output/\n")
cat("============================================================\n")
cat(" Figures:\n")
cat("   Figure_PCA.pdf\n")
cat("   Figure_FST_heatmap.pdf\n")
cat("   Figure_AMOVA.pdf\n")
cat("   Figure_Heterozygosity.pdf\n")
cat("   Figure_UPGMA_tree.pdf\n")
cat("   Figure_MSN.pdf\n")
cat("   Figure_DAPC.pdf\n")
cat("   Figure_PCoA.pdf\n")
cat(" Tables:\n")
cat("   PCA_eigenvalues.csv\n")
cat("   FST_matrix.txt\n")
cat("   FST_pvalues.txt\n")
cat("   Heterozygosity_summary.csv\n")
cat("   AMOVA_results.txt\n")
cat("============================================================\n")
cat(" Note: sNMF admixture analysis is run separately via LEA.\n")
cat("   See 03_admixture_LEA.R\n")
cat("============================================================\n")
