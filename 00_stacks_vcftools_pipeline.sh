#!/usr/bin/env bash
# =============================================================================
# 00_stacks_vcftools_pipeline.sh
# De novo assembly and SNP filtering pipeline for C. femoratus GBS data
#
# Species  : Cheirimedon femoratus (amphipod)
# Study    : Population genetics across South Shetland Islands
# Populations: DICF, KGCF, LICF, SICF
#
# Pipeline steps:
#   1. De novo mapping and SNP calling with Stacks (denovo_map.pl)
#   2. SNP filtering with VCFtools:
#      a. Missing data, MAF, and minimum depth filter
#      b. Remove individuals with >80% missing data
#      c. Minimum mean depth filter
#
# Usage:
#   bash 00_stacks_vcftools_pipeline.sh
#
# Requirements:
#   Stacks >= 2.0    (https://catchenlab.life.illinois.edu/stacks/)
#   VCFtools >= 0.1.15 (https://vcftools.github.io/)
#
# Input files:
#   samples/        — directory containing demultiplexed FASTQ files
#   Amphi.txt       — population map (sample_id  population)
#
# Output:
#   70whole/        — Stacks output directory
#   final.recode.vcf — final filtered VCF for downstream analyses
# =============================================================================

set -euo pipefail

# ── Parameters ────────────────────────────────────────────────────────────────
SAMPLES_DIR="samples"
POPMAP="Amphi.txt"
STACKS_OUTDIR="70whole"
MAX_MISSING=0.8       # minimum genotyping rate (--max-missing)
MAF=0.04              # minimum minor allele frequency
MIN_DP=3              # minimum per-individual depth per locus
MAX_HET=0.6           # maximum observed heterozygosity (Stacks)
MIND_THRESH=0.80      # remove individuals with >80% missing data
MIN_MEAN_DP=5         # minimum mean depth across all loci per individual

echo "============================================================"
echo " Stacks + VCFtools SNP filtering pipeline"
echo " Species: Cheirimedon femoratus"
echo "============================================================"

# ── Step 1: De novo assembly and SNP calling with Stacks ─────────────────────
# Parameters:
#   -M 2   : mismatches allowed between stacks within individuals
#   -n 3   : mismatches allowed between stacks between individuals
#   -p 0.8 : locus must be present in at least 80% of populations
#   --max-obs-het 0.6 : filter loci with observed heterozygosity > 0.6
#   --write-random-snp : retain one random SNP per locus (reduces linkage)
#   --vcf  : output VCF format

echo ""
echo "[1/4] Running Stacks denovo_map.pl..."

denovo_map.pl \
    --samples "$SAMPLES_DIR" \
    --paired \
    --popmap "$POPMAP" \
    -o "$STACKS_OUTDIR" \
    -M 2 \
    -n 3 \
    -p 0.8 \
    -X "populations: --max-obs-het $MAX_HET --write-random-snp --vcf"

echo "      Stacks complete. Output: $STACKS_OUTDIR/populations.snps.vcf"

# ── Step 2a: SNP-level filters ────────────────────────────────────────────────
# --max-missing 0.8 : retain SNPs present in >= 80% of individuals
# --maf 0.04        : minimum minor allele frequency
# --minDP 3         : minimum genotype depth

echo ""
echo "[2/4] Applying SNP-level filters (missing=$MAX_MISSING, MAF=$MAF, minDP=$MIN_DP)..."

vcftools \
    --vcf "$STACKS_OUTDIR/populations.snps.vcf" \
    --max-missing "$MAX_MISSING" \
    --maf "$MAF" \
    --minDP "$MIN_DP" \
    --recode \
    --recode-INFO-all \
    --out miss80maf4dp3

echo "      SNPs retained: $(grep -v '^#' miss80maf4dp3.recode.vcf | wc -l)"

# ── Step 2b: Identify and remove low-quality individuals ─────────────────────
# Calculate per-individual missingness, then remove individuals
# with more than 80% missing genotypes across retained loci

echo ""
echo "[3/4] Identifying individuals with >${MIND_THRESH}% missing data..."

vcftools \
    --vcf miss80maf4dp3.recode.vcf \
    --missing-indv \
    --out miss80maf4dp3

# Report missingness
echo "      Per-individual missing data summary:"
cat miss80maf4dp3.imiss

# Extract low-quality individuals
awk -v thresh="$MIND_THRESH" 'NR>1 && $5 > thresh {print $1}' \
    miss80maf4dp3.imiss > lowDP-80.indv

N_REMOVE=$(wc -l < lowDP-80.indv)
echo "      Individuals to remove (>${MIND_THRESH} missing): $N_REMOVE"

if [[ "$N_REMOVE" -gt 0 ]]; then
    echo "      Removing:"
    cat lowDP-80.indv
    vcftools \
        --vcf miss80maf4dp3.recode.vcf \
        --remove lowDP-80.indv \
        --recode \
        --recode-INFO-all \
        --out miss80maf4dp3INDV
else
    echo "      No individuals removed."
    cp miss80maf4dp3.recode.vcf miss80maf4dp3INDV.recode.vcf
fi

echo "      Individuals retained: $(bcftools query -l miss80maf4dp3INDV.recode.vcf | wc -l)"

# ── Step 2c: Minimum mean depth filter ───────────────────────────────────────
# Retain only loci where the mean sequencing depth per individual >= 5x

echo ""
echo "[4/4] Applying minimum mean depth filter (>=${MIN_MEAN_DP}x)..."

vcftools \
    --vcf miss80maf4dp3INDV.recode.vcf \
    --min-meanDP "$MIN_MEAN_DP" \
    --recode \
    --recode-INFO-all \
    --out final

# Rename to standard output name
mv final.recode.vcf final.recode.vcf

echo ""
echo "============================================================"
echo " Filtering complete."
echo "============================================================"
echo " Final VCF     : final.recode.vcf"
echo " Final SNPs    : $(grep -v '^#' final.recode.vcf | wc -l)"
echo " Final samples : $(bcftools query -l final.recode.vcf | wc -l)"
echo ""
echo " Next steps:"
echo "   bash 01_prepare_data.sh --vcf final.recode.vcf --popmap $POPMAP"
echo "   Rscript 02_population_genetics.R"
echo "   Rscript Ne_estimation_LDmethod.R"
echo "============================================================"
