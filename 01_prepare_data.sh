#!/usr/bin/env bash
# =============================================================================
# 01_prepare_data.sh
# Converts VCF to per-population PLINK PED/MAP files for Ne estimation
#
# Species  : Cheirimedon femoratus (amphipod)
# Study    : Population genetics across South Shetland Islands
#
# Input    : VCF file + population map (popmap)
# Output   : Per-population PED/MAP files in ne_data/populations/
#
# Usage    :
#   bash 01_prepare_data.sh --vcf final.recode.vcf --popmap Amphi.txt
#   bash 01_prepare_data.sh --vcf final.recode.vcf --popmap Amphi.txt \
#                           --outdir ne_data --maf 0.05 --geno 0.10 --mind 0.20
#
# Requirements:
#   PLINK v1.9  (https://www.cog-genomics.org/plink/)
#   module load PLINK/1.09b6.16-GCC-9.2.0  (NeSI)
#
# Popmap format (tab or space separated, no header):
#   sample_id   POPULATION
#   DICF_01     DICF
#   KGCF_01     KGCF
# =============================================================================

set -euo pipefail

# ── Default parameters ────────────────────────────────────────────────────────
VCF=""
POPMAP=""
OUTDIR="ne_data"
MAF=0.05       # minimum minor allele frequency
GENO=0.10      # maximum missing genotypes per SNP
MIND=0.20      # maximum missing data per individual
CHRSET=95      # max autosomes (PLINK 1.9 integer only — for non-human organisms)

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
  echo ""
  echo "Usage: $0 --vcf <file.vcf> --popmap <popmap.txt> [options]"
  echo ""
  echo "Required:"
  echo "  --vcf      Path to input VCF file"
  echo "  --popmap   Path to population map file (sample_id  population)"
  echo ""
  echo "Optional:"
  echo "  --outdir   Output directory (default: ne_data)"
  echo "  --maf      Minimum minor allele frequency (default: 0.05)"
  echo "  --geno     Maximum missing genotypes per SNP (default: 0.10)"
  echo "  --mind     Maximum missing data per individual (default: 0.20)"
  echo ""
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf)    VCF="$2";    shift 2 ;;
    --popmap) POPMAP="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --maf)    MAF="$2";    shift 2 ;;
    --geno)   GENO="$2";   shift 2 ;;
    --mind)   MIND="$2";   shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown argument: $1"; usage ;;
  esac
done

# ── Pre-flight checks ─────────────────────────────────────────────────────────
[[ -z "$VCF"    ]] && { echo "ERROR: --vcf is required";    usage; }
[[ -z "$POPMAP" ]] && { echo "ERROR: --popmap is required"; usage; }
[[ -f "$VCF"    ]] || { echo "ERROR: VCF not found: $VCF";       exit 1; }
[[ -f "$POPMAP" ]] || { echo "ERROR: popmap not found: $POPMAP"; exit 1; }
command -v plink &>/dev/null || {
  echo "ERROR: plink not in PATH"
  echo "       On NeSI run: module load PLINK/1.09b6.16-GCC-9.2.0"
  exit 1
}

# ── Setup output directories ──────────────────────────────────────────────────
rm -rf "$OUTDIR/plink_all"
mkdir -p "$OUTDIR/logs" "$OUTDIR/plink_all" "$OUTDIR/populations"

echo "============================================================"
echo " Data preparation for Ne estimation"
echo " Species: Cheirimedon femoratus"
echo "============================================================"
echo " VCF    : $VCF"
echo " popmap : $POPMAP"
echo " outdir : $OUTDIR"
echo " MAF=$MAF  GENO=$GENO  MIND=$MIND"
echo "------------------------------------------------------------"

# ── Step 1: VCF → PLINK binary ───────────────────────────────────────────────
echo ""
echo "[1/4] Converting VCF to PLINK binary format..."

plink \
  --vcf "$VCF" \
  --make-bed \
  --out "$OUTDIR/plink_all/all_samples" \
  --allow-extra-chr \
  --chr-set $CHRSET \
  --double-id \
  --set-missing-var-ids '@:#' \
  2>&1 | tee "$OUTDIR/logs/step1_vcf_convert.log"

[[ -f "$OUTDIR/plink_all/all_samples.bed" ]] || {
  echo "ERROR: VCF conversion failed. Check: $OUTDIR/logs/step1_vcf_convert.log"
  exit 1
}

RAW_SNPS=$(wc -l < "$OUTDIR/plink_all/all_samples.bim")
RAW_IND=$(wc -l  < "$OUTDIR/plink_all/all_samples.fam")
echo "  Variants loaded : $RAW_SNPS"
echo "  Individuals     : $RAW_IND"

# ── Step 2: Quality filtering ─────────────────────────────────────────────────
echo ""
echo "[2/4] Applying quality filters (MAF=$MAF, GENO=$GENO, MIND=$MIND)..."

plink \
  --bfile "$OUTDIR/plink_all/all_samples" \
  --maf "$MAF" \
  --geno "$GENO" \
  --mind "$MIND" \
  --make-bed \
  --out "$OUTDIR/plink_all/all_filtered" \
  --allow-extra-chr \
  --chr-set $CHRSET \
  2>&1 | tee "$OUTDIR/logs/step2_filter.log"

[[ -f "$OUTDIR/plink_all/all_filtered.bed" ]] || {
  echo "ERROR: Filtering failed. Check: $OUTDIR/logs/step2_filter.log"
  exit 1
}

FILT_SNPS=$(wc -l < "$OUTDIR/plink_all/all_filtered.bim")
FILT_IND=$(wc -l  < "$OUTDIR/plink_all/all_filtered.fam")
echo "  SNPs before filtering : $RAW_SNPS"
echo "  SNPs after  filtering : $FILT_SNPS"
echo "  Individuals retained  : $FILT_IND"

# ── Step 3: Verify popmap sample IDs match FAM file ──────────────────────────
echo ""
echo "[3/4] Checking popmap IDs against PLINK FAM file..."

FAM_IDS=$(cut -d' ' -f2 "$OUTDIR/plink_all/all_filtered.fam" | sort)
POPMAP_IDS=$(awk '{print $1}' "$POPMAP" | sort)
MATCHED=$(comm -12 <(echo "$FAM_IDS") <(echo "$POPMAP_IDS") | wc -l)
ONLY_FAM=$(comm -23 <(echo "$FAM_IDS") <(echo "$POPMAP_IDS") | wc -l)
ONLY_POP=$(comm -13 <(echo "$FAM_IDS") <(echo "$POPMAP_IDS") | wc -l)

echo "  IDs in both FAM and popmap : $MATCHED"
echo "  IDs in FAM but not popmap  : $ONLY_FAM  (will be skipped)"
echo "  IDs in popmap but not FAM  : $ONLY_POP  (check for typos if unexpected)"

if [[ "$MATCHED" -eq 0 ]]; then
  echo ""
  echo "ERROR: No sample IDs matched between FAM and popmap!"
  echo "First 5 FAM IDs:"
  cut -d' ' -f2 "$OUTDIR/plink_all/all_filtered.fam" | head -5
  echo "First 5 popmap IDs:"
  awk '{print $1}' "$POPMAP" | head -5
  echo "IDs must match exactly (same spelling, case, and separators)."
  exit 1
fi

# ── Step 4: Split by population and produce PED/MAP files ────────────────────
echo ""
echo "[4/4] Splitting into per-population PED/MAP files..."

POPULATIONS=$(awk '{print $2}' "$POPMAP" | sort -u)
echo "  Populations found: $(echo $POPULATIONS | tr '\n' ' ')"

for POP in $POPULATIONS; do
  echo ""
  echo "  ── $POP ────────────────────────────────"

  POPDIR="$OUTDIR/populations/$POP"
  mkdir -p "$POPDIR"

  # Build keep file: FID IID (both set to sample_id via --double-id)
  awk -v pop="$POP" '$2 == pop {print $1, $1}' "$POPMAP" > "$POPDIR/keep.txt"
  N_SAMPLES=$(wc -l < "$POPDIR/keep.txt")
  echo "  Samples in popmap : $N_SAMPLES"

  if [[ "$N_SAMPLES" -eq 0 ]]; then
    echo "  WARNING: No samples found for $POP — skipping."
    continue
  fi

  if [[ "$N_SAMPLES" -lt 10 ]]; then
    echo "  WARNING: Only $N_SAMPLES samples — Ne estimates may be unreliable."
  fi

  # Extract population → binary BED
  plink \
    --bfile "$OUTDIR/plink_all/all_filtered" \
    --keep "$POPDIR/keep.txt" \
    --maf 0.01 \
    --make-bed \
    --out "$POPDIR/${POP}_binary" \
    --allow-extra-chr \
    --chr-set $CHRSET \
    2>&1 | tee "$OUTDIR/logs/step4_${POP}_binary.log"

  # Binary → PED/MAP (text format required by Ne_estimation_LDmethod.R)
  plink \
    --bfile "$POPDIR/${POP}_binary" \
    --recode \
    --out "$POPDIR/${POP}" \
    --allow-extra-chr \
    --chr-set $CHRSET \
    2>&1 | tee "$OUTDIR/logs/step4_${POP}_ped.log"

  if [[ -f "$POPDIR/${POP}.ped" ]]; then
    N_IND=$(wc -l < "$POPDIR/${POP}.ped")
    N_SNP=$(wc -l < "$POPDIR/${POP}.map")
    echo "  Individuals : $N_IND"
    echo "  SNPs        : $N_SNP"
    echo "  Output      : $POPDIR/${POP}.ped / .map"
  else
    echo "  WARNING: PED file not created. Check: $OUTDIR/logs/step4_${POP}_ped.log"
  fi
done

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " Summary"
echo "============================================================"
printf "  %-10s  %-12s  %-10s\n" "Population" "Individuals" "SNPs"
printf "  %-10s  %-12s  %-10s\n" "----------" "-----------" "----"
for POP in $POPULATIONS; do
  POPDIR="$OUTDIR/populations/$POP"
  if [[ -f "$POPDIR/${POP}.ped" ]]; then
    printf "  %-10s  %-12s  %-10s\n" \
      "$POP" \
      "$(wc -l < $POPDIR/${POP}.ped)" \
      "$(wc -l < $POPDIR/${POP}.map)"
  else
    printf "  %-10s  %-12s  %-10s\n" "$POP" "FAILED" "—"
  fi
done
echo ""
echo " Next step:"
echo "   Rscript Ne_estimation_LDmethod.R --datadir $OUTDIR/populations \\"
echo "                                     --outdir  ne_results"
echo "============================================================"
