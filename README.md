# Amphipod-population-genetics-across-South-Shetland-Islands
We investigated the population genomics of the common benthic brooding Antarctic amphipod Cheirimedon femoratus,  from 87 individuals collected across four sites in the South Shetland Islands: Deception Island (DICF), King George Island (KGCF), Livingston Island (LICF), and Snow Island (SICF)

########## In Bash

# SNP Calling Script(Amphipod)

module load Stacks
denovo_map.pl --samples samples/ --paired --popmap Amphipod.txt -o test1 -M 2 -m 3 -n 3 -r 0.8 -X "populations: --max-obs-het 0.6 --write-random-snp --vcf"

# Filter based on missing data and minor allele frequency
module load VCFtools
vcftools --vcf populations.snps.vcf --min-meanDP 4 --max-meanDP 60 --recode --recode-INFO-all
mv out.recode.vcf depth_filtered.vcf
vcftools --vcf depth_filtered.vcf --max-missing 0.4 --recode --recode-INFO-all
mv out.recode.vcf depth_maxmiss40.vcf



############## Visualization in R

setwd("/Users/sepva552/Desktop/GBS/Amphipod)

# Load required libraries
library(vcfR)
library(adegenet)
library(StAMPP)
library(ggplot2)
library(corrplot)
library(LEA)
library(pophelper)
library(reshape2)
library(diveRsity)
library(devtools)
library(vegan)
library(SNPRelate)
library(gdsfmt)
library(ape)
library(phangorn)
library(parallel)
library(dartR)
library(pegas)

##############################
snp_vcf2 = read.vcfR("./final.recode.vcf")
pop.data2 = read.table("./Amphi.txt", header = F)
gl.snp2 <- vcfR2genlight(snp_vcf2)
pop(gl.snp2) <- rep(pop.data2$V2)
snp.pca2 <- glPca(gl.snp2, nf = 10)
snp.pca.scores2 <- as.data.frame(snp.pca2$scores)
snp.pca.scores2$pop <- pop(gl.snp2)
write.table(snp.pca.scores2, "./Amphi_PCA.txt", sep = "\t")
eig.val<-snp.pca2$eig
eig.val
nInd(gl.snp2)
length(pop.data2$V2) 
#percentages of variance for each PC
eig.perc <- 100*snp.pca2$eig/sum(snp.pca2$eig)
eig.perc
eigen<-data.frame(eig.val,eig.perc)
eigen

##PCA plotting
data2 = read.delim("Amphi_adegenetPCA.txt") #I have manually added population information to this file prior to loading it
mycol = c("#f1c039", "#f37d21", "#51692d", "#56ba32", "#87CEFA")
ggplot(data2, aes(x=PC1, y=PC2, color=pop)) +
  geom_point(size = 3) + 
  scale_color_manual(values=mycol) +
  theme_classic()+
  xlab("PC1 (13.89%)") +
  ylab("PC2 (2.21%)")


##calculating Fst between populations
Qfly_Fst <- stamppFst(gl.snp2, nboots = 100, percent = 95, nclusters = 4)
Fst <- Qfly_Fst$Fsts
pFst <- Qfly_Fst$Pvalues
write.table(Fst, "Fst1.txt", sep="\t")
write.table(pFst, "Fst1_pvalue.txt", sep="\t")
##creating heatmap
# Melt the correlation matrix
library(ggplot2)
library(reshape2)

# Define the FST matrix
fst_matrix <- matrix(
  c(NA, NA, NA, NA,
    0.00528, NA, NA, NA,
    0.01219, 0.01716, NA, NA,
    0.00539, 0.00727, 0.01109, NA),
  nrow = 4, byrow = TRUE,
  dimnames = list(c("DICF", "KGCF", "LICF", "SICF"),
                  c("DICF", "KGCF", "LICF", "SICF"))
)

# Melt the matrix for ggplot
melted_fst <- melt(fst_matrix, na.rm = TRUE)

# Reverse the order of Var2 (for left-side labels)
melted_fst$Var2 <- factor(melted_fst$Var2, levels = rev(levels(melted_fst$Var2)))

# Plot with y-axis labels on the left
ggplot(melted_fst, aes(Var2, Var1, fill = value)) +  # Note: Var2 and Var1 swapped
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 5)), size = 3.5) +
  scale_fill_gradient(low = "blue", high = "yellow", na.value = "white",
                      limits = c(0, 0.02), name = "FST") +
  labs(title = "Pairwise FST between South Shetland amphipod populations",
       x = "", y = "") +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1))  # Align y-axis labels to the left


library(LEA)
library(pophelper)

##creating input files
vcf2geno(input.file = "final.recode.vcf", output.file = "Qff.geno")

##snmf clustering
projectalpha = NULL
projectalpha = snmf("Qff.geno", K = 1:10, repetitions = 50, entropy = T, CPU = 8, project = "new")

# plot cross-entropy criterion for all runs in the snmf project
pdf(file = "./cross_ent_alphadefualt.pdf")
plot(projectalpha, col = "maroon4", pch = 19, cex = 1.2)
dev.off()

best2 = which.min(cross.entropy(projectalpha, K = 2))
best2
best3 = which.min(cross.entropy(projectalpha, K = 3))
best3
best4 = which.min(cross.entropy(projectalpha, K = 4))
best4
best5 = which.min(cross.entropy(projectalpha, K = 5))
best5
best6 = which.min(cross.entropy(projectalpha, K = 6))
best6
best7 = which.min(cross.entropy(projectalpha, K = 7))
best7
best8 = which.min(cross.entropy(projectalpha, K = 8))
best8
best9 = which.min(cross.entropy(projectalpha, K = 9))
best9
best10 = which.min(cross.entropy(projectalpha, K = 10))
best10


##creating admixture plots. For this, you need to first create a new folder (All_Qfiles) and move the Q files with "best" entropies from the LEA runs into it. 
sfiles <- list.files(path=("./All_Qfiles"), full.names=T)
slist <- readQ(files=sfiles)
clustercol <- c("#51692d", "#f1c039", "#f37d21", "#56ba32", "#a1d99b",
                "#9ecae1", "#fc9272", "#dd1c77", "#756bb1", "#636363")
plotQ(qlist=slist[3],imgtype = "pdf",
      height = 1.5, clustercol = c("#51692d","#56ba32","#f1c039"), dpi = 1200, exportpath = "./")
plotQ(qlist=slist[4],imgtype = "pdf",
      height = 1.5, clustercol = c("#51692d","#56ba32","#f1c039","#f37d21"), dpi = 1200, exportpath = "./")
plotQ(qlist=slist[5],imgtype = "pdf",
      height = 1.5, clustercol = c("#f37d21","#51692d","#f1c039","#56ba32","#a63838"), dpi = 1200, exportpath = "./")
plotQ(qlist=slist[6],imgtype = "pdf",
      height = 1.5, clustercol = c("#a63838","#f1c039","#ecbcab","#56ba32","#51692d","#f37d21"), dpi = 1200, exportpath = "./")
plotQ(qlist=slist[7],imgtype = "pdf",
      height = 1.5, clustercol = c("#51692d","#a63838","#f1c039","#ecbcab","#caf291","#56ba32","#f37d21"), dpi = 1200, exportpath = "./")
plotQ(qlist=slist[8],imgtype = "pdf",
      height = 1.5, clustercol = c("#3d87db","#a63838","#51692d","#f1c039","#ecbcab","#caf291","#56ba32","#f37d21"), dpi = 1200, exportpath = "./")
plotQ(qlist=slist[9],imgtype = "pdf",
      height = 1.5, clustercol = c("#f1c039","#a63838","#3d87db","#ecbcab","#caf291","#56ba32","#0000FF","#51692d","#f37d21"), dpi = 1200, exportpath = "./")
plotQ(qlist=slist[1],imgtype = "pdf",
      height = 1.5, clustercol = c("#dd00ff","#51692d","#caf291","#3d87db","#ecbcab","#a63838","#56ba32","#0000FF","#f1c039","#f37d21"), dpi = 1200, exportpath = "./") 






library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)



# Extract cross-entropy values
Kvals <- 1:10
cv_errors <- sapply(Kvals, function(k) {
  min(cross.entropy(projectalpha, K = k))
})

# Create a data frame for plotting
cv_df <- data.frame(
  K = Kvals,
  CrossEntropy = cv_errors
)

# Publication-style plot
cv_plot <- ggplot(cv_df, aes(x = K, y = CrossEntropy)) +
  geom_line(color = "black", size = 0.8) +
  geom_point(color = "maroon4", size = 3) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    x = "Number of ancestral populations (K)",
    y = "Minimum cross-entropy",
    title = "Cross-validation error plot for sNMF (LEA)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )

# Save high-resolution version
ggsave("CrossEntropy_CV_Plot_LEA.tiff", cv_plot, dpi = 600, width = 12, height = 9, units = "cm", compression = "lzw")
ggsave("CrossEntropy_CV_Plot_LEA.pdf", cv_plot, width = 12, height = 9, units = "cm")

# Print the plot
print(cv_plot)


########## Isolation By Distance

# FST matrix (symmetric)
fst_matrix <- matrix(c(
  NA, 0.0052776956, 0.0121943661, 0.0053894006,
  NA, NA, 0.0171592160, 0.0072681676,
  NA, NA, NA, 0.0110918219,
  NA, NA, NA, NA
), nrow = 4, byrow = TRUE,
dimnames = list(c("DICF", "KGCF", "LICF", "SICF"),
                c("DICF", "KGCF", "LICF", "SICF")))
fst_matrix[lower.tri(fst_matrix)] <- t(fst_matrix)[lower.tri(fst_matrix)]

# GPS data
gps_data <- data.frame(
  site = c("KGCF", "LICF", "SICF", "DICF"),
  longitude = c(-58.780489, -60.621406, -61.209067, -60.709267),
  latitude  = c(-62.220081, -62.643944, -62.748483, -62.968181)
)

# Compute distance matrices
geo_dist <- as.dist(distm(gps_data[, c("longitude", "latitude")], fun = distGeo) / 1000)
fst_dist <- as.dist(fst_matrix)
mantel_result <- mantel(geo_dist, fst_dist, permutations = 999)

# Prepare data for plotting
pairwise_data <- melt(fst_matrix, na.rm = TRUE) %>%
  filter(as.character(Var1) < as.character(Var2)) %>%
  rename(site1 = Var1, site2 = Var2, FST = value) %>%
  mutate(
    distance = distGeo(
      cbind(gps_data$longitude[match(site1, gps_data$site)],
            gps_data$latitude[match(site1, gps_data$site)]),
      cbind(gps_data$longitude[match(site2, gps_data$site)],
            gps_data$latitude[match(site2, gps_data$site)])
    ) / 1000,
    pair_label = paste(site1, site2, sep = "–")
  )

# Fit linear model
lm_model <- lm(FST ~ distance, data = pairwise_data)
lm_summary <- summary(lm_model)

# Create plot
mantel_plot <- ggplot(pairwise_data, aes(x = distance, y = FST)) +
  geom_point(size = 3.5, color = "black") +
  geom_text_repel(aes(label = pair_label), size = 3, segment.color = "grey50") +
  geom_smooth(method = "lm", color = "black", fill = "lightgray", se = TRUE, size = 0.8, linetype = "solid") +
  labs(
    x = "Geographic distance (km)",
    y = expression(F[ST]),
    caption = sprintf("Mantel test: r = %.3f, p = %.3f\nLinear model: R² = %.3f, slope = %.2e, p = %.3f",
                      mantel_result$statistic, mantel_result$signif,
                      lm_summary$r.squared, lm_summary$coefficients[2, 1],
                      lm_summary$coefficients[2, 4])
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.caption = element_text(hjust = 0, size = 10)
  )

# Save high-res versions
ggsave("Fig_IBD_Mantel_MolEcol.tiff", plot = mantel_plot, dpi = 600, width = 16, height = 12, units = "cm", compression = "lzw")
ggsave("Fig_IBD_Mantel_MolEcol.pdf", plot = mantel_plot, width = 16, height = 12, units = "cm")

# Print statistical summary
cat("=== Mantel and Linear Model Results ===\n")
cat(sprintf("Mantel Test: r = %.3f, p = %.3f\n", mantel_result$statistic, mantel_result$signif))
cat(sprintf("Linear Regression: R² = %.3f, slope = %.2e, p = %.3f\n", lm_summary$r.squared, lm_summary$coefficients[2, 1], lm_summary$coefficients[2, 4]))


##############AMOVA
library(vcfR)       # Read VCF
library(adegenet)   # Convert to genind
library(poppr)      # AMOVA
library(dartR)      # Quality checks (optional)
install.packages("dartR")
# Read VCF
vcf <- read.vcfR("final.recode.vcf")

# Read population assignments
pop_data <- read.table("Amphi.txt", header = FALSE, col.names = c("Sample", "Population"))

# Convert VCF to genind
genind_data <- vcfR2genind(vcf)

# Assign populations (order must match VCF samples!)
pop(genind_data) <- pop_data$Population[match(indNames(genind_data), pop_data$Sample)]

# Verify
table(pop(genind_data))

genind_data <- genind_data[pop(genind_data) != "Pop1", ]
# Basic AMOVA (no hierarchy)
amova_result <- poppr.amova(genind_data, ~Population, within = FALSE)

# Test significance with 1000 permutations
amova_signif <- randtest(amova_result, nrepet = 1000)
plot(amova_signif)  # Check p-value distribution
amova_signif

# Scale genotypes
genind_scaled <- scaleGen(genind_data, NA.method = "mean")

# Run PCA
pca <- dudi.pca(genind_scaled, cent = FALSE, scale = FALSE, scannf = FALSE)

# Plot
library(ggplot2)
ggplot(data.frame(pca$li, Population = pop(genind_data)), 
       aes(x = Axis1, y = Axis2, color = Population)) +
  geom_point(size = 3) +
  stat_ellipse() +  # Optional: Add ellipses
  theme_minimal()
dapc_result <- dapc(genind_data, n.pca = 10, n.da = 2)  # Adjust n.pca as needed
scatter(dapc_result, col = rainbow(3))




# Install required packages if not installed(He and Ho)


# Load libraries
library(adegenet)
library(vcfR)
library(hierfstat)
install.packages("hierfstat")

library(poppr)
library(ggplot2)
library(reshape2)  # For data reshaping

# 1. Read VCF and assign populations (same as before)
vcf <- read.vcfR("final.recode.vcf")
genind_obj <- vcfR2genind(vcf)
poptx <- read.table("Amphi.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(poptx) <- c("SampleID", "Population")
pop(genind_obj) <- poptx$Population[match(indNames(genind_obj), poptx$SampleID)]
hierfstat_data <- genind2hierfstat(genind_obj)

# 2. Calculate Ho and He
hierfstat_data <- genind2hierfstat(genind_obj)

basic_stats <- basic.stats(hierfstat_data, diploid = TRUE)
Ho <- basic_stats$Ho
He <- basic_stats$Hs



# 3. Plot Ho vs. He
heterozygosity_data <- data.frame(
  Population = colnames(Ho),
  Ho = colMeans(Ho, na.rm = TRUE),
  He = colMeans(He, na.rm = TRUE)
)

# Reshape for ggplot
heterozygosity_melted <- melt(heterozygosity_data, id.vars = "Population")

# Plot
ggplot(heterozygosity_melted, aes(x = Population, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  labs(title = "Observed vs. Expected Heterozygosity",
       y = "Heterozygosity", x = "Population", fill = "Metric") +
  scale_fill_manual(values = c("Ho" = "#4E79A7", "He" = "#F28E2B")) +
  theme_minimal() +
  theme(legend.position = "top")

######## Folded SFS
# Load required libraries
library(vcfR)
library(adegenet)
library(pegas)

# Define file paths
vcf_file <- "final.recode.vcf"      # Your input VCF file
popmap_file <- "Amphi.txt"          # Two-column file: sampleID, population

# Read the VCF
vcf <- read.vcfR(vcf_file)

# Convert to genind object
gen <- vcfR2genind(vcf)

# Read the population map
popmap <- read.table(popmap_file, header = FALSE, col.names = c("sample", "pop"))

# Match sample names and assign populations
inds <- indNames(gen)
pop_assignments <- popmap$pop[match(inds, popmap$sample)]
gen@pop <- as.factor(pop_assignments)

# Set layout: 2 rows × 2 columns
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))  # adjust margins if needed

# Function to compute and plot folded SFS
compute_folded_sfs <- function(subgen, pop_name) {
  ac <- tab(subgen, NA.method = "mean")
  
  if (all(is.na(ac)) || nrow(ac) == 0) {
    message(paste("Skipping", pop_name, "- no valid genotypes"))
    return(NULL)
  }
  
  maf <- apply(ac, 2, function(x) {
    if (all(is.na(x))) return(NA)
    af <- mean(x, na.rm = TRUE) / 2
    min(af, 1 - af)
  })
  
  maf <- maf[!is.na(maf) & maf >= 0 & maf <= 0.5]
  
  message(paste(pop_name, "- valid SNPs:", length(maf)))
  
  if (length(maf) < 5) {
    message(paste("Skipping", pop_name, "- too few SNPs"))
    return(NULL)
  }
  
  hist(maf,
       breaks = seq(0, 0.5, by = 0.05),
       main = paste("Folded SFS -", pop_name),
       xlab = "Minor Allele Frequency",
       ylab = "SNP Count",
       col = "skyblue",
       border = "gray")
}

# Loop through each population and plot SFS
unique_pops <- unique(pop(gen))
for (p in unique_pops) {
  subgen <- gen[pop(gen) == p]
  compute_folded_sfs(subgen, as.character(p))
}

# Optional: reset layout
par(mfrow = c(1, 1))





######################### Tajima D

##Plotting Sliding windows analysis - TajimaD
library(ggplot2)
library(dplyr)
library(forcats)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(ggpubr)
library(rstatix)


##########  in bash
##module load VCFtools
## Create a Population File (for example Pop1)
### vcftools --vcf final.recode.vcf --keep Pop1.txt --TajimaD 10000 --out pop1_tajimaD
###vcftools --vcf final.recode.vcf --keep Pop2.txt --TajimaD 10000 --out pop2_tajimaD
###vcftools --vcf final.recode.vcf --keep Pop3.txt --TajimaD 10000 --out pop3_tajimaD
##### for Nucleotide Diversity (π).   vcftools --vcf final.recode.vcf --keep Pop1.txt --window-pi 10000 --out Pop1_pi



# 1. Load and clean data
LICF <- read.table("tajima_LICF.Tajima.D", header = TRUE)
LICF <- LICF[!is.na(LICF$TajimaD) & LICF$TajimaD != "NaN", ]
LICF$pop <- "LICF"

DICF <- read.table("tajima_DICF.Tajima.D", header = TRUE)
DICF <- DICF[!is.na(DICF$TajimaD) & DICF$TajimaD != "NaN", ]
DICF$pop <- "DICF"

SICF <- read.table("tajima_SICF.Tajima.D", header = TRUE)
SICF <- SICF[!is.na(SICF$TajimaD) & SICF$TajimaD != "NaN", ]
SICF$pop <- "SICF"

KGCF <- read.table("tajima_KGCF.Tajima.D", header = TRUE)
KGCF <- KGCF[!is.na(KGCF$TajimaD) & KGCF$TajimaD != "NaN", ]
KGCF$pop <- "KGCF"

# 2. Combine into one dataframe
all_data <- rbind(LICF, DICF, SICF, KGCF)
all_data$pop <- factor(all_data$pop, levels = c("KGCF", "DICF", "LICF", "SICF"))

library(ggplot2)
library(viridis)

ggplot(all_data, aes(x = pop, y = TajimaD, fill = pop)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(
    title = "Tajima’s D distribution across populations",
    x = "Population",
    y = "Tajima's D"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )



##### Statistics
KGCFTajimaD <- read.table("tajima_KGCF.Tajima.D", header = TRUE)

# Clean data: remove NA or NaN values
KGCFTajimaD <- KGCFTajimaD[!is.na(KGCFTajimaD$TajimaD) &KGCFTajimaD$TajimaD != "NaN", ]

# Convert to numeric just in case
KGCFTajimaD$TajimaD <- as.numeric(KGCFTajimaD$TajimaD)

# 1. Summary statistics
summary_stats <- summary(KGCFTajimaD$TajimaD)
print(summary_stats)

# 2. Median
med <- median(KGCFTajimaD$TajimaD)
cat("Median Tajima's D:", med, "\n")

# 3. Mean
mean_val <- mean(KGCFTajimaD$TajimaD)
cat("Mean Tajima's D:", mean_val, "\n")

# 4. Standard Deviation
sd_val <- sd(KGCFTajimaD$TajimaD)
cat("Standard Deviation:", sd_val, "\n")

# 5. 5th and 95th Percentiles
quant_vals <- quantile(KGCFTajimaD$TajimaD, probs = c(0.05, 0.95))
cat("5th Percentile:", quant_vals[1], "\n")
cat("95th Percentile:", quant_vals[2], "\n")

# 6. Min and Max (optional)
cat("Min:", min(KGCFTajimaD$TajimaD), "\n")
cat("Max:", max(KGCFTajimaD$TajimaD), "\n")



###### Nucleotide diversity

#########IN bash
vcftools --vcf final.recode.vcf --keep pop1.txt --window-pi 10000 --window-pi-step 10000 --out pop1_windowed_pi
##Repeat for pop2 and pop3 and pop4 

KGCF_pi <- fread("KGCF.windowed.pi", header = TRUE)
LICF_pi <- fread("LICF.windowed.pi", header = TRUE)
SICF_pi <- fread("SICF.windowed.pi", header = TRUE)
DICF_pi <- fread("DICF.windowed.pi", header = TRUE)

# Prepare each dataset: remove NA, keep only 'PI' column and label population
KGCF_pi_clean <- KGCF_pi %>%
  filter(!is.na(PI)) %>%
  select(PI) %>%
  mutate(pop = "KGCF")

LICF_pi_clean <- LICF_pi %>%
  filter(!is.na(PI)) %>%
  select(PI) %>%
  mutate(pop = "LICF")

SICF_pi_clean <- SICF_pi %>%
  filter(!is.na(PI)) %>%
  select(PI) %>%
  mutate(pop = "SICF")

DICF_pi_clean <- DICF_pi %>%
  filter(!is.na(PI)) %>%
  select(PI) %>%
  mutate(pop = "DICF")

# Combine all into one dataframe
pi_all <- rbind(KGCF_pi_clean, LICF_pi_clean, SICF_pi_clean, DICF_pi_clean)

# Confirm data
table(pi_all$pop)

# Violin + boxplot
ggplot(pi_all, aes(x = pop, y = PI, fill = pop)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Nucleotide Diversity (π) Across Amphipod Populations",
    x = "Population",
    y = expression(pi)
  )

# Kruskal-Wallis test
kruskal_result <- kruskal.test(PI ~ pop, data = pi_all)
print(kruskal_result)

# Pairwise Wilcoxon post-hoc test
pairwise_result <- pairwise.wilcox.test(pi_all$PI, pi_all$pop, p.adjust.method = "bonferroni")

print(pairwise_result)



######TREEMIX
#####treemIX
####IN BASH:
  populations -V ./final.recode.vcf -O ./ -M ./Amphi.txt --treemix
gzip final.recode.p.treemix
module load Miniconda3
for i in {1..10}; do treemix -i final.recode.p.treemix.gz  -m ${i} -bootstrap -k 1000 -o BMSB_${i} > treemix_${i}_log ; done
#####in R
source("plotting_funcs.R") 
plot_tree("BMSB_8")
plot_tree("BMSB_3")
plot_tree("BMSB_10")



