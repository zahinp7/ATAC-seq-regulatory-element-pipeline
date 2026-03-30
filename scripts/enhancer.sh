#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 04:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=9

module load bedtools

# Set paths
HUMAN_GTF="/path/to/your/base/directory/HumanGenomeInfo/gencode.v47.annotation.gff3.gz"
MOUSE_GTF="/path/to/your/base/directory/MouseGenomeInfo/gencode.vM10.annotation.gff3.gz"
OUTPUT_DIR="/path/to/your/base/directory/classified_regions"

# Set peak file paths (Input for classification)
HUMAN_LIVER_PEAKS="/path/to/your/base/directory/HumanAtac/Liver/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"
HUMAN_PANCREAS_PEAKS="/path/to/your/base/directory/HumanAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"
MOUSE_LIVER_PEAKS="/path/to/your/base/directory/MouseAtac/Liver/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"
MOUSE_PANCREAS_PEAKS="/path/to/your/base/directory/MouseAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"

# Create base output directory
mkdir -p $OUTPUT_DIR

# ===== STEP 3: Classify Enhancers vs Promoters & Generate Subsets =====
# Note: This combined step replaces the previous separate classification, comparison, and statistics sections.
echo "Starting Classification and Subset Generation..."
# Create directories for the results of classification and comparisons
CLASSIFIED_DIR="$OUTPUT_DIR/classified_peaks"
COMP_DIR="$OUTPUT_DIR/region_comparisons" # Changed name to reflect promoters+enhancers
mkdir -p "$CLASSIFIED_DIR/human" "$CLASSIFIED_DIR/mouse"
mkdir -p "$COMP_DIR/tissue_comparison" "$COMP_DIR/species_comparison"

# --- PART 3.1: Define Promoter Regions ---
HUMAN_PROMOTERS_DEF_BED="$CLASSIFIED_DIR/human_promoters_definition.bed" # File defining promoter coords
MOUSE_PROMOTERS_DEF_BED="$CLASSIFIED_DIR/mouse_promoters_definition.bed" # File defining promoter coords

echo "Generating promoter definition BED files..."
zcat "$HUMAN_GTF" | awk 'BEGIN{OFS="\t"} $3=="transcript" {if ($7=="+"){start=$4-2000;if(start<0)start=0;print $1,start,$4+200,"promoter",".",$7}else{start=$5-200;if(start<0)start=0;print $1,start,$5+2000,"promoter",".",$7}}' | grep -v "^chrM" > "$HUMAN_PROMOTERS_DEF_BED"
zcat "$MOUSE_GTF" | awk 'BEGIN{OFS="\t"} $3=="transcript" {if ($7=="+"){start=$4-2000;if(start<0)start=0;print $1,start,$4+200,"promoter",".",$7}else{start=$5-200;if(start<0)start=0;print $1,start,$5+2000,"promoter",".",$7}}' | grep -v "^chrM" > "$MOUSE_PROMOTERS_DEF_BED"

# --- PART 3.2: Classify Original Peaks ---
echo "Classifying original peaks..."

# Define output file paths for initial classification (peaks overlapping/not overlapping definitions)
HUMAN_LIVER_PROMOTERS_BED="$CLASSIFIED_DIR/human/human_liver_promoters.bed"
HUMAN_LIVER_ENHANCERS_BED="$CLASSIFIED_DIR/human/human_liver_enhancers.bed"
HUMAN_PANCREAS_PROMOTERS_BED="$CLASSIFIED_DIR/human/human_pancreas_promoters.bed"
HUMAN_PANCREAS_ENHANCERS_BED="$CLASSIFIED_DIR/human/human_pancreas_enhancers.bed"
MOUSE_LIVER_PROMOTERS_BED="$CLASSIFIED_DIR/mouse/mouse_liver_promoters.bed"
MOUSE_LIVER_ENHANCERS_BED="$CLASSIFIED_DIR/mouse/mouse_liver_enhancers.bed"
MOUSE_PANCREAS_PROMOTERS_BED="$CLASSIFIED_DIR/mouse/mouse_pancreas_promoters.bed"
MOUSE_PANCREAS_ENHANCERS_BED="$CLASSIFIED_DIR/mouse/mouse_pancreas_enhancers.bed"

# Function to classify peaks (Promoter vs Enhancer)
classify_peaks() {
  local input_peak_gz=$1 # Expecting the original gzipped narrowPeak
  local promoter_def_bed=$2 # The promoter *definition* file
  local output_promoter_peaks_bed=$3 # Peaks overlapping definitions
  local output_enhancer_peaks_bed=$4 # Peaks NOT overlapping definitions
  local tmp_peak_bed="/tmp/$(basename "$input_peak_gz" .narrowPeak.gz)_$$_peaks.bed"

  echo "  Classifying $(basename "$input_peak_gz") ..."
  zcat "$input_peak_gz" | cut -f 1-3 > "$tmp_peak_bed"
  bedtools intersect -a "$tmp_peak_bed" -b "$promoter_def_bed" -u > "$output_promoter_peaks_bed"
  bedtools intersect -a "$tmp_peak_bed" -b "$promoter_def_bed" -v > "$output_enhancer_peaks_bed"
  rm "$tmp_peak_bed"
  echo "  Finished classifying $(basename "$input_peak_gz")"
}

# Run classification
classify_peaks "$HUMAN_LIVER_PEAKS" "$HUMAN_PROMOTERS_DEF_BED" "$HUMAN_LIVER_PROMOTERS_BED" "$HUMAN_LIVER_ENHANCERS_BED"
classify_peaks "$HUMAN_PANCREAS_PEAKS" "$HUMAN_PROMOTERS_DEF_BED" "$HUMAN_PANCREAS_PROMOTERS_BED" "$HUMAN_PANCREAS_ENHANCERS_BED"
classify_peaks "$MOUSE_LIVER_PEAKS" "$MOUSE_PROMOTERS_DEF_BED" "$MOUSE_LIVER_PROMOTERS_BED" "$MOUSE_LIVER_ENHANCERS_BED"
classify_peaks "$MOUSE_PANCREAS_PEAKS" "$MOUSE_PROMOTERS_DEF_BED" "$MOUSE_PANCREAS_PROMOTERS_BED" "$MOUSE_PANCREAS_ENHANCERS_BED"

# --- PART 3.3: Tissue Comparison of Enhancers ---
echo "Performing tissue comparison on enhancer sets..."

# Define output file paths
HUMAN_ENHANCERS_TISSUE_SHARED="$COMP_DIR/tissue_comparison/human_enhancers_tissue_shared.bed"
HUMAN_LIVER_ENHANCERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/human_liver_enhancers_tissue_specific.bed"
HUMAN_PANCREAS_ENHANCERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/human_pancreas_enhancers_tissue_specific.bed"
MOUSE_ENHANCERS_TISSUE_SHARED="$COMP_DIR/tissue_comparison/mouse_enhancers_tissue_shared.bed"
MOUSE_LIVER_ENHANCERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/mouse_liver_enhancers_tissue_specific.bed"
MOUSE_PANCREAS_ENHANCERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/mouse_pancreas_enhancers_tissue_specific.bed"

# Human Enhancer Tissue Comparison
bedtools intersect -a "$HUMAN_LIVER_ENHANCERS_BED" -b "$HUMAN_PANCREAS_ENHANCERS_BED" -u > "$HUMAN_ENHANCERS_TISSUE_SHARED"
bedtools intersect -a "$HUMAN_LIVER_ENHANCERS_BED" -b "$HUMAN_PANCREAS_ENHANCERS_BED" -v > "$HUMAN_LIVER_ENHANCERS_TISSUE_SPECIFIC"
bedtools intersect -a "$HUMAN_PANCREAS_ENHANCERS_BED" -b "$HUMAN_LIVER_ENHANCERS_BED" -v > "$HUMAN_PANCREAS_ENHANCERS_TISSUE_SPECIFIC"

# Mouse Enhancer Tissue Comparison
bedtools intersect -a "$MOUSE_LIVER_ENHANCERS_BED" -b "$MOUSE_PANCREAS_ENHANCERS_BED" -u > "$MOUSE_ENHANCERS_TISSUE_SHARED"
bedtools intersect -a "$MOUSE_LIVER_ENHANCERS_BED" -b "$MOUSE_PANCREAS_ENHANCERS_BED" -v > "$MOUSE_LIVER_ENHANCERS_TISSUE_SPECIFIC"
bedtools intersect -a "$MOUSE_PANCREAS_ENHANCERS_BED" -b "$MOUSE_LIVER_ENHANCERS_BED" -v > "$MOUSE_PANCREAS_ENHANCERS_TISSUE_SPECIFIC"

# --- PART 3.4: Tissue Comparison of Promoters ---
echo "Performing tissue comparison on promoter sets..."

# Define output file paths
HUMAN_PROMOTERS_TISSUE_SHARED="$COMP_DIR/tissue_comparison/human_promoters_tissue_shared.bed"
HUMAN_LIVER_PROMOTERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/human_liver_promoters_tissue_specific.bed"
HUMAN_PANCREAS_PROMOTERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/human_pancreas_promoters_tissue_specific.bed"
MOUSE_PROMOTERS_TISSUE_SHARED="$COMP_DIR/tissue_comparison/mouse_promoters_tissue_shared.bed"
MOUSE_LIVER_PROMOTERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/mouse_liver_promoters_tissue_specific.bed"
MOUSE_PANCREAS_PROMOTERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/mouse_pancreas_promoters_tissue_specific.bed"

# Human Promoter Tissue Comparison
bedtools intersect -a "$HUMAN_LIVER_PROMOTERS_BED" -b "$HUMAN_PANCREAS_PROMOTERS_BED" -u > "$HUMAN_PROMOTERS_TISSUE_SHARED"
bedtools intersect -a "$HUMAN_LIVER_PROMOTERS_BED" -b "$HUMAN_PANCREAS_PROMOTERS_BED" -v > "$HUMAN_LIVER_PROMOTERS_TISSUE_SPECIFIC"
bedtools intersect -a "$HUMAN_PANCREAS_PROMOTERS_BED" -b "$HUMAN_LIVER_PROMOTERS_BED" -v > "$HUMAN_PANCREAS_PROMOTERS_TISSUE_SPECIFIC"

# Mouse Promoter Tissue Comparison
bedtools intersect -a "$MOUSE_LIVER_PROMOTERS_BED" -b "$MOUSE_PANCREAS_PROMOTERS_BED" -u > "$MOUSE_PROMOTERS_TISSUE_SHARED"
bedtools intersect -a "$MOUSE_LIVER_PROMOTERS_BED" -b "$MOUSE_PANCREAS_PROMOTERS_BED" -v > "$MOUSE_LIVER_PROMOTERS_TISSUE_SPECIFIC"
bedtools intersect -a "$MOUSE_PANCREAS_PROMOTERS_BED" -b "$MOUSE_LIVER_PROMOTERS_BED" -v > "$MOUSE_PANCREAS_PROMOTERS_TISSUE_SPECIFIC"

# --- PART 3.5: Species Comparison (Enhancers vs Mapped Peaks) ---
echo "Performing species comparison (enhancers vs mapped peaks)..."
# IMPORTANT: Define the location of the *processed* mapped *peak* files (BED format) generated by HALPER/Step 2 of the pipeline
MAPPED_PEAKS_PROCESSED_DIR="/path/to/your/base/directory/results/cross_species/processed_files" # Example path, adjust!

MAPPED_HL_TO_MM_PEAKS="$MAPPED_PEAKS_PROCESSED_DIR/human_liver_to_mouse.bed"
MAPPED_HP_TO_MM_PEAKS="$MAPPED_PEAKS_PROCESSED_DIR/human_pancreas_to_mouse.bed"
MAPPED_ML_TO_HG_PEAKS="$MAPPED_PEAKS_PROCESSED_DIR/mouse_liver_to_human.bed"
MAPPED_MP_TO_HG_PEAKS="$MAPPED_PEAKS_PROCESSED_DIR/mouse_pancreas_to_human.bed"

# Check if mapped peak files exist before attempting comparison
SKIP_SPECIES_STATS=0 # Default to 0 (don't skip)
if [[ ! -f "$MAPPED_ML_TO_HG_PEAKS" || ! -f "$MAPPED_MP_TO_HG_PEAKS" || ! -f "$MAPPED_HL_TO_MM_PEAKS" || ! -f "$MAPPED_HP_TO_MM_PEAKS" ]]; then
  echo "Warning: One or more processed mapped peak files not found in $MAPPED_PEAKS_PROCESSED_DIR. Skipping species comparison steps."
  SKIP_SPECIES_STATS=1
fi

if [[ "$SKIP_SPECIES_STATS" -eq 0 ]]; then
  # Define output file paths for Enhancer comparison
  HUMAN_LIVER_ENH_SHARED_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_liver_enhancers_shared_with_mapped_mouse_peaks_human_coords.bed"
  MOUSE_LIVER_ENH_SHARED_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_liver_enhancers_shared_with_mapped_human_peaks_mouse_coords.bed"
  HUMAN_LIVER_ENH_SPECIFIC_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_liver_enhancers_specific_vs_mapped_mouse_peaks_human_coords.bed"
  MOUSE_LIVER_ENH_SPECIFIC_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_liver_enhancers_specific_vs_mapped_human_peaks_mouse_coords.bed"
  HUMAN_PANCREAS_ENH_SHARED_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_pancreas_enhancers_shared_with_mapped_mouse_peaks_human_coords.bed"
  MOUSE_PANCREAS_ENH_SHARED_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_pancreas_enhancers_shared_with_mapped_human_peaks_mouse_coords.bed"
  HUMAN_PANCREAS_ENH_SPECIFIC_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_pancreas_enhancers_specific_vs_mapped_mouse_peaks_human_coords.bed"
  MOUSE_PANCREAS_ENH_SPECIFIC_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_pancreas_enhancers_specific_vs_mapped_human_peaks_mouse_coords.bed"

  # Liver Species Comparison (Enhancers vs Mapped Peaks)
  bedtools intersect -a "$HUMAN_LIVER_ENHANCERS_BED" -b "$MAPPED_ML_TO_HG_PEAKS" -u > "$HUMAN_LIVER_ENH_SHARED_MAPPED_MOUSE_HG_COORDS"
  bedtools intersect -a "$MOUSE_LIVER_ENHANCERS_BED" -b "$MAPPED_HL_TO_MM_PEAKS" -u > "$MOUSE_LIVER_ENH_SHARED_MAPPED_HUMAN_MM_COORDS"
  bedtools intersect -a "$HUMAN_LIVER_ENHANCERS_BED" -b "$MAPPED_ML_TO_HG_PEAKS" -v > "$HUMAN_LIVER_ENH_SPECIFIC_MAPPED_MOUSE_HG_COORDS"
  bedtools intersect -a "$MOUSE_LIVER_ENHANCERS_BED" -b "$MAPPED_HL_TO_MM_PEAKS" -v > "$MOUSE_LIVER_ENH_SPECIFIC_MAPPED_HUMAN_MM_COORDS"

  # Pancreas Species Comparison (Enhancers vs Mapped Peaks)
  bedtools intersect -a "$HUMAN_PANCREAS_ENHANCERS_BED" -b "$MAPPED_MP_TO_HG_PEAKS" -u > "$HUMAN_PANCREAS_ENH_SHARED_MAPPED_MOUSE_HG_COORDS"
  bedtools intersect -a "$MOUSE_PANCREAS_ENHANCERS_BED" -b "$MAPPED_HP_TO_MM_PEAKS" -u > "$MOUSE_PANCREAS_ENH_SHARED_MAPPED_HUMAN_MM_COORDS"
  bedtools intersect -a "$HUMAN_PANCREAS_ENHANCERS_BED" -b "$MAPPED_MP_TO_HG_PEAKS" -v > "$HUMAN_PANCREAS_ENH_SPECIFIC_MAPPED_MOUSE_HG_COORDS"
  bedtools intersect -a "$MOUSE_PANCREAS_ENHANCERS_BED" -b "$MAPPED_HP_TO_MM_PEAKS" -v > "$MOUSE_PANCREAS_ENH_SPECIFIC_MAPPED_HUMAN_MM_COORDS"
fi

# --- PART 3.6: Species Comparison (Promoters vs Mapped Peaks) ---
echo "Performing species comparison (promoters vs mapped peaks)..."

if [[ "$SKIP_SPECIES_STATS" -eq 0 ]]; then
  # Define output file paths for Promoter comparison
  HUMAN_LIVER_PROM_SHARED_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_liver_promoters_shared_with_mapped_mouse_peaks_human_coords.bed"
  MOUSE_LIVER_PROM_SHARED_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_liver_promoters_shared_with_mapped_human_peaks_mouse_coords.bed"
  HUMAN_LIVER_PROM_SPECIFIC_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_liver_promoters_specific_vs_mapped_mouse_peaks_human_coords.bed"
  MOUSE_LIVER_PROM_SPECIFIC_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_liver_promoters_specific_vs_mapped_human_peaks_mouse_coords.bed"
  HUMAN_PANCREAS_PROM_SHARED_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_pancreas_promoters_shared_with_mapped_mouse_peaks_human_coords.bed"
  MOUSE_PANCREAS_PROM_SHARED_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_pancreas_promoters_shared_with_mapped_human_peaks_mouse_coords.bed"
  HUMAN_PANCREAS_PROM_SPECIFIC_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_pancreas_promoters_specific_vs_mapped_mouse_peaks_human_coords.bed"
  MOUSE_PANCREAS_PROM_SPECIFIC_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_pancreas_promoters_specific_vs_mapped_human_peaks_mouse_coords.bed"

  # Liver Species Comparison (Promoters vs Mapped Peaks)
  bedtools intersect -a "$HUMAN_LIVER_PROMOTERS_BED" -b "$MAPPED_ML_TO_HG_PEAKS" -u > "$HUMAN_LIVER_PROM_SHARED_MAPPED_MOUSE_HG_COORDS"
  bedtools intersect -a "$MOUSE_LIVER_PROMOTERS_BED" -b "$MAPPED_HL_TO_MM_PEAKS" -u > "$MOUSE_LIVER_PROM_SHARED_MAPPED_HUMAN_MM_COORDS"
  bedtools intersect -a "$HUMAN_LIVER_PROMOTERS_BED" -b "$MAPPED_ML_TO_HG_PEAKS" -v > "$HUMAN_LIVER_PROM_SPECIFIC_MAPPED_MOUSE_HG_COORDS"
  bedtools intersect -a "$MOUSE_LIVER_PROMOTERS_BED" -b "$MAPPED_HL_TO_MM_PEAKS" -v > "$MOUSE_LIVER_PROM_SPECIFIC_MAPPED_HUMAN_MM_COORDS"

  # Pancreas Species Comparison (Promoters vs Mapped Peaks)
  bedtools intersect -a "$HUMAN_PANCREAS_PROMOTERS_BED" -b "$MAPPED_MP_TO_HG_PEAKS" -u > "$HUMAN_PANCREAS_PROM_SHARED_MAPPED_MOUSE_HG_COORDS"
  bedtools intersect -a "$MOUSE_PANCREAS_PROMOTERS_BED" -b "$MAPPED_HP_TO_MM_PEAKS" -u > "$MOUSE_PANCREAS_PROM_SHARED_MAPPED_HUMAN_MM_COORDS"
  bedtools intersect -a "$HUMAN_PANCREAS_PROMOTERS_BED" -b "$MAPPED_MP_TO_HG_PEAKS" -v > "$HUMAN_PANCREAS_PROM_SPECIFIC_MAPPED_MOUSE_HG_COORDS"
  bedtools intersect -a "$MOUSE_PANCREAS_PROMOTERS_BED" -b "$MAPPED_HP_TO_MM_PEAKS" -v > "$MOUSE_PANCREAS_PROM_SPECIFIC_MAPPED_HUMAN_MM_COORDS"
fi

# --- PART 3.7: Calculate Statistics & Percentages ---
echo "Calculating statistics and percentages..."
STATS_FILE="$OUTPUT_DIR/classification_and_comparison_summary.txt"

# Function to safely calculate percentages using bc
calc_pct() {
    local num=$1
    local denom=$2
    if [[ "$denom" -eq 0 || -z "$denom" ]]; then
        echo "0.00" # Return 0 if denominator is zero or empty
    else
        # Use bc for floating point arithmetic, scale=2 for 2 decimal places
        echo "scale=2; ($num / $denom) * 100" | bc
    fi
}

# Get all the counts first
# Initial Classification
HLP_C=$(wc -l < "$HUMAN_LIVER_PROMOTERS_BED")
HLE_C=$(wc -l < "$HUMAN_LIVER_ENHANCERS_BED")
HPP_C=$(wc -l < "$HUMAN_PANCREAS_PROMOTERS_BED")
HPE_C=$(wc -l < "$HUMAN_PANCREAS_ENHANCERS_BED")
MLP_C=$(wc -l < "$MOUSE_LIVER_PROMOTERS_BED")
MLE_C=$(wc -l < "$MOUSE_LIVER_ENHANCERS_BED")
MPP_C=$(wc -l < "$MOUSE_PANCREAS_PROMOTERS_BED")
MPE_C=$(wc -l < "$MOUSE_PANCREAS_ENHANCERS_BED")
# Enhancer Tissue Comparison
HETS_C=$(wc -l < "$HUMAN_ENHANCERS_TISSUE_SHARED")
HLETS_C=$(wc -l < "$HUMAN_LIVER_ENHANCERS_TISSUE_SPECIFIC")
HPETS_C=$(wc -l < "$HUMAN_PANCREAS_ENHANCERS_TISSUE_SPECIFIC")
METS_C=$(wc -l < "$MOUSE_ENHANCERS_TISSUE_SHARED")
MLETS_C=$(wc -l < "$MOUSE_LIVER_ENHANCERS_TISSUE_SPECIFIC")
MPETS_C=$(wc -l < "$MOUSE_PANCREAS_ENHANCERS_TISSUE_SPECIFIC")
# Promoter Tissue Comparison
HPTS_C=$(wc -l < "$HUMAN_PROMOTERS_TISSUE_SHARED")
HPLTS_C=$(wc -l < "$HUMAN_LIVER_PROMOTERS_TISSUE_SPECIFIC")
HPPTS_C=$(wc -l < "$HUMAN_PANCREAS_PROMOTERS_TISSUE_SPECIFIC")
MPTS_C=$(wc -l < "$MOUSE_PROMOTERS_TISSUE_SHARED")
MPLTS_C=$(wc -l < "$MOUSE_LIVER_PROMOTERS_TISSUE_SPECIFIC")
MPPTS_C=$(wc -l < "$MOUSE_PANCREAS_PROMOTERS_TISSUE_SPECIFIC")

# Species Comparison Counts (initialize to 0, update if not skipped)
HLES_C=0; MLES_C=0; HLESPEC_C=0; MLESPEC_C=0; HPES_C=0; MPES_C=0; HPESPEC_C=0; MPESPEC_C=0
HPLS_C=0; MPLS_C=0; HPLSPEC_C=0; MPLSPEC_C=0; HPPS_C=0; MPPS_C=0; HPPSPEC_C=0; MPPSPEC_C=0
if [[ "$SKIP_SPECIES_STATS" -eq 0 ]]; then
  HLES_C=$(wc -l < "$HUMAN_LIVER_ENH_SHARED_MAPPED_MOUSE_HG_COORDS")
  MLES_C=$(wc -l < "$MOUSE_LIVER_ENH_SHARED_MAPPED_HUMAN_MM_COORDS")
  HLESPEC_C=$(wc -l < "$HUMAN_LIVER_ENH_SPECIFIC_MAPPED_MOUSE_HG_COORDS")
  MLESPEC_C=$(wc -l < "$MOUSE_LIVER_ENH_SPECIFIC_MAPPED_HUMAN_MM_COORDS")
  HPES_C=$(wc -l < "$HUMAN_PANCREAS_ENH_SHARED_MAPPED_MOUSE_HG_COORDS")
  MPES_C=$(wc -l < "$MOUSE_PANCREAS_ENH_SHARED_MAPPED_HUMAN_MM_COORDS")
  HPESPEC_C=$(wc -l < "$HUMAN_PANCREAS_ENH_SPECIFIC_MAPPED_MOUSE_HG_COORDS")
  MPESPEC_C=$(wc -l < "$MOUSE_PANCREAS_ENH_SPECIFIC_MAPPED_HUMAN_MM_COORDS")
  HPLS_C=$(wc -l < "$HUMAN_LIVER_PROM_SHARED_MAPPED_MOUSE_HG_COORDS")
  MPLS_C=$(wc -l < "$MOUSE_LIVER_PROM_SHARED_MAPPED_HUMAN_MM_COORDS")
  HPLSPEC_C=$(wc -l < "$HUMAN_LIVER_PROM_SPECIFIC_MAPPED_MOUSE_HG_COORDS")
  MPLSPEC_C=$(wc -l < "$MOUSE_LIVER_PROM_SPECIFIC_MAPPED_HUMAN_MM_COORDS")
  HPPS_C=$(wc -l < "$HUMAN_PANCREAS_PROM_SHARED_MAPPED_MOUSE_HG_COORDS")
  MPPS_C=$(wc -l < "$MOUSE_PANCREAS_PROM_SHARED_MAPPED_HUMAN_MM_COORDS")
  HPPSPEC_C=$(wc -l < "$HUMAN_PANCREAS_PROM_SPECIFIC_MAPPED_MOUSE_HG_COORDS")
  MPPSPEC_C=$(wc -l < "$MOUSE_PANCREAS_PROM_SPECIFIC_MAPPED_HUMAN_MM_COORDS")
fi

# Calculate Percentages
# Tissue Sharing
HL_ENH_SHARED_PCT=$(calc_pct $HETS_C $HLE_C)
HP_ENH_SHARED_PCT=$(calc_pct $HETS_C $HPE_C)
HL_PROM_SHARED_PCT=$(calc_pct $HPTS_C $HLP_C)
HP_PROM_SHARED_PCT=$(calc_pct $HPTS_C $HPP_C)
ML_ENH_SHARED_PCT=$(calc_pct $METS_C $MLE_C)
MP_ENH_SHARED_PCT=$(calc_pct $METS_C $MPE_C)
ML_PROM_SHARED_PCT=$(calc_pct $MPTS_C $MLP_C)
MP_PROM_SHARED_PCT=$(calc_pct $MPTS_C $MPP_C)
# Species Sharing (vs Mapped Peaks)
HL_ENH_SHARED_SPECIES_PCT=$(calc_pct $HLES_C $HLE_C)
HP_ENH_SHARED_SPECIES_PCT=$(calc_pct $HPES_C $HPE_C)
HL_PROM_SHARED_SPECIES_PCT=$(calc_pct $HPLS_C $HLP_C)
HP_PROM_SHARED_SPECIES_PCT=$(calc_pct $HPPS_C $HPP_C)
ML_ENH_SHARED_SPECIES_PCT=$(calc_pct $MLES_C $MLE_C)
MP_ENH_SHARED_SPECIES_PCT=$(calc_pct $MPES_C $MPE_C)
ML_PROM_SHARED_SPECIES_PCT=$(calc_pct $MPLS_C $MLP_C)
MP_PROM_SHARED_SPECIES_PCT=$(calc_pct $MPPS_C $MPP_C)


# Write statistics file
echo "===== Classification and Comparison Summary =====" > "$STATS_FILE"
echo "" >> "$STATS_FILE"
echo "--- Initial Classification Counts ---" >> "$STATS_FILE"
echo "Human liver promoters: $HLP_C" >> "$STATS_FILE"
echo "Human liver enhancers: $HLE_C" >> "$STATS_FILE"
echo "Human pancreas promoters: $HPP_C" >> "$STATS_FILE"
echo "Human pancreas enhancers: $HPE_C" >> "$STATS_FILE"
echo "Mouse liver promoters: $MLP_C" >> "$STATS_FILE"
echo "Mouse liver enhancers: $MLE_C" >> "$STATS_FILE"
echo "Mouse pancreas promoters: $MPP_C" >> "$STATS_FILE"
echo "Mouse pancreas enhancers: $MPE_C" >> "$STATS_FILE"
echo "" >> "$STATS_FILE"

echo "--- Question 1: Tissue Sharing Comparison (Counts & Percentages) ---" >> "$STATS_FILE"
echo "** Human **" >> "$STATS_FILE"
echo "  Shared Enhancers (Liver/Pancreas): $HETS_C" >> "$STATS_FILE"
echo "  Shared Promoters (Liver/Pancreas): $HPTS_C" >> "$STATS_FILE"
echo "  Liver Enhancers Shared: ${HL_ENH_SHARED_PCT}% ($HETS_C/$HLE_C)" >> "$STATS_FILE"
echo "  Liver Promoters Shared: ${HL_PROM_SHARED_PCT}% ($HPTS_C/$HLP_C)" >> "$STATS_FILE"
echo "  Pancreas Enhancers Shared: ${HP_ENH_SHARED_PCT}% ($HETS_C/$HPE_C)" >> "$STATS_FILE"
echo "  Pancreas Promoters Shared: ${HP_PROM_SHARED_PCT}% ($HPTS_C/$HPP_C)" >> "$STATS_FILE"
echo "** Mouse **" >> "$STATS_FILE"
echo "  Shared Enhancers (Liver/Pancreas): $METS_C" >> "$STATS_FILE"
echo "  Shared Promoters (Liver/Pancreas): $MPTS_C" >> "$STATS_FILE"
echo "  Liver Enhancers Shared: ${ML_ENH_SHARED_PCT}% ($METS_C/$MLE_C)" >> "$STATS_FILE"
echo "  Liver Promoters Shared: ${ML_PROM_SHARED_PCT}% ($MPTS_C/$MLP_C)" >> "$STATS_FILE"
echo "  Pancreas Enhancers Shared: ${MP_ENH_SHARED_PCT}% ($METS_C/$MPE_C)" >> "$STATS_FILE"
echo "  Pancreas Promoters Shared: ${MP_PROM_SHARED_PCT}% ($MPTS_C/$MPP_C)" >> "$STATS_FILE"
echo "" >> "$STATS_FILE"

echo "--- Question 2: Species Sharing Comparison vs Mapped Peaks (Counts & Percentages) ---" >> "$STATS_FILE"
if [[ "$SKIP_SPECIES_STATS" -eq 1 ]]; then
  echo "Skipped due to missing mapped peak files." >> "$STATS_FILE"
else
  echo "Note: Sharing means overlap with orthologs of *any* peak from the other species." >> "$STATS_FILE"
  echo "** Human Liver vs Mouse Mapped Peaks **" >> "$STATS_FILE"
  echo "  Shared Enhancers (Human Coords): $HLES_C" >> "$STATS_FILE"
  echo "  Shared Promoters (Human Coords): $HPLS_C" >> "$STATS_FILE"
  echo "  Human Liver Enhancers Shared: ${HL_ENH_SHARED_SPECIES_PCT}% ($HLES_C/$HLE_C)" >> "$STATS_FILE"
  echo "  Human Liver Promoters Shared: ${HL_PROM_SHARED_SPECIES_PCT}% ($HPLS_C/$HLP_C)" >> "$STATS_FILE"
  echo "** Human Pancreas vs Mouse Mapped Peaks **" >> "$STATS_FILE"
  echo "  Shared Enhancers (Human Coords): $HPES_C" >> "$STATS_FILE"
  echo "  Shared Promoters (Human Coords): $HPPS_C" >> "$STATS_FILE"
  echo "  Human Pancreas Enhancers Shared: ${HP_ENH_SHARED_SPECIES_PCT}% ($HPES_C/$HPE_C)" >> "$STATS_FILE"
  echo "  Human Pancreas Promoters Shared: ${HP_PROM_SHARED_SPECIES_PCT}% ($HPPS_C/$HPP_C)" >> "$STATS_FILE"
  echo "** Mouse Liver vs Human Mapped Peaks **" >> "$STATS_FILE"
  echo "  Shared Enhancers (Mouse Coords): $MLES_C" >> "$STATS_FILE"
  echo "  Shared Promoters (Mouse Coords): $MPLS_C" >> "$STATS_FILE"
  echo "  Mouse Liver Enhancers Shared: ${ML_ENH_SHARED_SPECIES_PCT}% ($MLES_C/$MLE_C)" >> "$STATS_FILE"
  echo "  Mouse Liver Promoters Shared: ${ML_PROM_SHARED_SPECIES_PCT}% ($MPLS_C/$MLP_C)" >> "$STATS_FILE"
  echo "** Mouse Pancreas vs Human Mapped Peaks **" >> "$STATS_FILE"
  echo "  Shared Enhancers (Mouse Coords): $MPES_C" >> "$STATS_FILE"
  echo "  Shared Promoters (Mouse Coords): $MPPS_C" >> "$STATS_FILE"
  echo "  Mouse Pancreas Enhancers Shared: ${MP_ENH_SHARED_SPECIES_PCT}% ($MPES_C/$MPE_C)" >> "$STATS_FILE"
  echo "  Mouse Pancreas Promoters Shared: ${MP_PROM_SHARED_SPECIES_PCT}% ($MPPS_C/$MPP_C)" >> "$STATS_FILE"
fi
echo "" >> "$STATS_FILE"

echo "Classification and comparisons complete. Summary in $STATS_FILE"