#!/bin/bash
# pro_enh.sh - Classify peaks into promoter/enhancer, compare subsets, calculate stats.
#
# Usage: ./pro_enh.sh <output_dir> <human_gtf> <mouse_gtf> <human_liver_peaks_gz> <human_pancreas_peaks_gz> <mouse_liver_peaks_gz> <mouse_pancreas_peaks_gz> <mapped_hl_to_mm_bed> <mapped_hp_to_mm_bed> <mapped_ml_to_hg_bed> <mapped_mp_to_hg_bed>
#
# Inputs:
#   Arg 1:  <output_dir>              - Base directory for pipeline results.
#   Arg 2:  <human_gtf>               - Gzipped GFF3: Human gene annotation.
#   Arg 3:  <mouse_gtf>               - Gzipped GFF3: Mouse gene annotation.
#   Arg 4:  <human_liver_peaks_gz>    - Gzipped narrowPeak: Original Human Liver peaks.
#   Arg 5:  <human_pancreas_peaks_gz> - Gzipped narrowPeak: Original Human Pancreas peaks.
#   Arg 6:  <mouse_liver_peaks_gz>    - Gzipped narrowPeak: Original Mouse Liver peaks.
#   Arg 7:  <mouse_pancreas_peaks_gz> - Gzipped narrowPeak: Original Mouse Pancreas peaks.
#   Arg 8:  <mapped_hl_to_mm_bed>     - BED: Mapped Human Liver peaks (Mouse coords, from Step 1).
#   Arg 9:  <mapped_hp_to_mm_bed>     - BED: Mapped Human Pancreas peaks (Mouse coords, from Step 1).
#   Arg 10: <mapped_ml_to_hg_bed>     - BED: Mapped Mouse Liver peaks (Human coords, from Step 1).
#   Arg 11: <mapped_mp_to_hg_bed>     - BED: Mapped Mouse Pancreas peaks (Human coords, from Step 1).
#
# Outputs:
#   Directory: <output_dir>/classification_and_comparisons/classified_peaks/
#   File:      .../human_promoters_definition.bed               - BED: Defined Human promoter regions (from GTF).
#   File:      .../mouse_promoters_definition.bed               - BED: Defined Mouse promoter regions (from GTF).
#   File:      .../human/human_liver_promoters.bed              - BED: Original Human liver peaks overlapping promoters.
#   File:      .../human/human_liver_enhancers.bed              - BED: Original Human liver peaks NOT overlapping promoters.
#   File:      .../human/human_pancreas_promoters.bed           - BED: Original Human pancreas peaks overlapping promoters.
#   File:      .../human/human_pancreas_enhancers.bed           - BED: Original Human pancreas peaks NOT overlapping promoters.
#   File:      .../mouse/mouse_liver_promoters.bed              - BED: Original Mouse liver peaks overlapping promoters.
#   File:      .../mouse/mouse_liver_enhancers.bed              - BED: Original Mouse liver peaks NOT overlapping promoters.
#   File:      .../mouse/mouse_pancreas_promoters.bed           - BED: Original Mouse pancreas peaks overlapping promoters.
#   File:      .../mouse/mouse_pancreas_enhancers.bed           - BED: Original Mouse pancreas peaks NOT overlapping promoters.
#
#   Directory: <output_dir>/classification_and_comparisons/region_comparisons/tissue_comparison/
#   File:      .../human_enhancers_tissue_shared.bed           - BED: Human enhancers shared between liver and pancreas.
#   File:      .../human_liver_enhancers_tissue_specific.bed     - BED: Human liver-specific enhancers.
#   File:      .../human_pancreas_enhancers_tissue_specific.bed  - BED: Human pancreas-specific enhancers.
#   File:      .../mouse_enhancers_tissue_shared.bed           - BED: Mouse enhancers shared between liver and pancreas.
#   File:      .../mouse_liver_enhancers_tissue_specific.bed     - BED: Mouse liver-specific enhancers.
#   File:      .../mouse_pancreas_enhancers_tissue_specific.bed  - BED: Mouse pancreas-specific enhancers.
#   File:      .../human_promoters_tissue_shared.bed           - BED: Human promoters shared between liver and pancreas.
#   File:      .../human_liver_promoters_tissue_specific.bed     - BED: Human liver-specific promoters.
#   File:      .../human_pancreas_promoters_tissue_specific.bed  - BED: Human pancreas-specific promoters.
#   File:      .../mouse_promoters_tissue_shared.bed           - BED: Mouse promoters shared between liver and pancreas.
#   File:      .../mouse_liver_promoters_tissue_specific.bed     - BED: Mouse liver-specific promoters.
#   File:      .../mouse_pancreas_promoters_tissue_specific.bed  - BED: Mouse pancreas-specific promoters.
#
#   Directory: <output_dir>/classification_and_comparisons/region_comparisons/species_comparison/
#   File:      .../human_liver_enhancers_shared_with_mapped_mouse_peaks_human_coords.bed    - BED: Human liver enhancers overlapping mapped mouse peaks.
#   File:      .../mouse_liver_enhancers_shared_with_mapped_human_peaks_mouse_coords.bed    - BED: Mouse liver enhancers overlapping mapped human peaks.
#   File:      .../human_liver_enhancers_specific_vs_mapped_mouse_peaks_human_coords.bed    - BED: Human liver enhancers NOT overlapping mapped mouse peaks.
#   File:      .../mouse_liver_enhancers_specific_vs_mapped_human_peaks_mouse_coords.bed    - BED: Mouse liver enhancers NOT overlapping mapped human peaks.
#   File:      .../human_pancreas_enhancers_shared_with_mapped_mouse_peaks_human_coords.bed - BED: Human pancreas enhancers overlapping mapped mouse peaks.
#   File:      .../mouse_pancreas_enhancers_shared_with_mapped_human_peaks_mouse_coords.bed - BED: Mouse pancreas enhancers overlapping mapped human peaks.
#   File:      .../human_pancreas_enhancers_specific_vs_mapped_mouse_peaks_human_coords.bed - BED: Human pancreas enhancers NOT overlapping mapped mouse peaks.
#   File:      .../mouse_pancreas_enhancers_specific_vs_mapped_human_peaks_mouse_coords.bed - BED: Mouse pancreas enhancers NOT overlapping mapped human peaks.
#   File:      .../human_liver_promoters_shared_with_mapped_mouse_peaks_human_coords.bed    - BED: Human liver promoters overlapping mapped mouse peaks.
#   File:      .../mouse_liver_promoters_shared_with_mapped_human_peaks_mouse_coords.bed    - BED: Mouse liver promoters overlapping mapped human peaks.
#   File:      .../human_liver_promoters_specific_vs_mapped_mouse_peaks_human_coords.bed    - BED: Human liver promoters NOT overlapping mapped mouse peaks.
#   File:      .../mouse_liver_promoters_specific_vs_mapped_human_peaks_mouse_coords.bed    - BED: Mouse liver promoters NOT overlapping mapped human peaks.
#   File:      .../human_pancreas_promoters_shared_with_mapped_mouse_peaks_human_coords.bed - BED: Human pancreas promoters overlapping mapped mouse peaks.
#   File:      .../mouse_pancreas_promoters_shared_with_mapped_human_peaks_mouse_coords.bed - BED: Mouse pancreas promoters overlapping mapped human peaks.
#   File:      .../human_pancreas_promoters_specific_vs_mapped_mouse_peaks_human_coords.bed - BED: Human pancreas promoters NOT overlapping mapped mouse peaks.
#   File:      .../mouse_pancreas_promoters_specific_vs_mapped_human_peaks_mouse_coords.bed - BED: Mouse pancreas promoters NOT overlapping mapped human peaks.
#
#   File: <output_dir>/classification_and_comparisons/classification_and_comparison_summary.txt - Text summary of classification and comparison statistics.

set -e

# Parse command-line arguments
if [ $# -lt 11 ]; then
    echo "Usage: $0 <output_dir> <human_gtf> <mouse_gtf> <human_liver_peaks_gz> <human_pancreas_peaks_gz> <mouse_liver_peaks_gz> <mouse_pancreas_peaks_gz> <mapped_hl_to_mm_bed> <mapped_hp_to_mm_bed> <mapped_ml_to_hg_bed> <mapped_mp_to_hg_bed>"
    exit 1
fi

OUTPUT_DIR=$1
HUMAN_GTF=$2
MOUSE_GTF=$3
HUMAN_LIVER_PEAKS=$4 # Original narrowPeak.gz path
HUMAN_PANCREAS_PEAKS=$5 # Original narrowPeak.gz path
MOUSE_LIVER_PEAKS=$6 # Original narrowPeak.gz path
MOUSE_PANCREAS_PEAKS=$7 # Original narrowPeak.gz path
MAPPED_HL_TO_MM_PEAKS=$8 # Mapped BED path from Step 1
MAPPED_HP_TO_MM_PEAKS=$9 # Mapped BED path from Step 1
MAPPED_ML_TO_HG_PEAKS=${10} # Mapped BED path from Step 1
MAPPED_MP_TO_HG_PEAKS=${11} # Mapped BED path from Step 1

# Define output directories based on main output directory
ANALYSIS_OUTPUT_DIR="$OUTPUT_DIR/classification_and_comparisons"
MAPPED_DIR="$OUTPUT_DIR/mapped_peaks" # Define where mapped peaks are expected (for SKIP_SPECIES_STATS check)

# Create necessary directories upfront
mkdir -p "$ANALYSIS_OUTPUT_DIR/classified_peaks/human" "$ANALYSIS_OUTPUT_DIR/classified_peaks/mouse"
mkdir -p "$ANALYSIS_OUTPUT_DIR/region_comparisons/tissue_comparison" "$ANALYSIS_OUTPUT_DIR/region_comparisons/species_comparison"

# ===== STEP 4: Classify Peaks, Compare Subsets & Calculate Stats =====
# This step generates the Promoter/Enhancer classifications and comparisons needed for the final questions and MEME-ChIP
echo "===== STEP 4: Classifying Peaks, Comparing Subsets & Calculating Stats ====="

# Define directories (already created above)
CLASSIFIED_DIR="$ANALYSIS_OUTPUT_DIR/classified_peaks"
COMP_DIR="$ANALYSIS_OUTPUT_DIR/region_comparisons"

# --- PART 4.1: Define Promoter Regions ---
HUMAN_PROMOTERS_DEF_BED="$CLASSIFIED_DIR/human_promoters_definition.bed"
MOUSE_PROMOTERS_DEF_BED="$CLASSIFIED_DIR/mouse_promoters_definition.bed"
echo "  Generating promoter definition BED files (if needed)..."
if [[ ! -s "$HUMAN_PROMOTERS_DEF_BED" ]]; then zcat "$HUMAN_GTF" | awk 'BEGIN{OFS="\t"} $3=="transcript" {if ($7=="+"){start=$4-2000;if(start<0)start=0;print $1,start,$4+200,"promoter",".",$7}else{start=$5-200;if(start<0)start=0;print $1,start,$5+2000,"promoter",".",$7}}' | grep -v "^chrM" > "$HUMAN_PROMOTERS_DEF_BED"; else echo "    Human promoter definition file already exists."; fi
if [[ ! -s "$MOUSE_PROMOTERS_DEF_BED" ]]; then zcat "$MOUSE_GTF" | awk 'BEGIN{OFS="\t"} $3=="transcript" {if ($7=="+"){start=$4-2000;if(start<0)start=0;print $1,start,$4+200,"promoter",".",$7}else{start=$5-200;if(start<0)start=0;print $1,start,$5+2000,"promoter",".",$7}}' | grep -v "^chrM" > "$MOUSE_PROMOTERS_DEF_BED"; else echo "    Mouse promoter definition file already exists."; fi

# --- PART 4.2: Classify Original Peaks ---
echo "  Classifying original peaks into Promoters/Enhancers..."
HUMAN_LIVER_PROMOTERS_BED="$CLASSIFIED_DIR/human/human_liver_promoters.bed"
HUMAN_LIVER_ENHANCERS_BED="$CLASSIFIED_DIR/human/human_liver_enhancers.bed"
HUMAN_PANCREAS_PROMOTERS_BED="$CLASSIFIED_DIR/human/human_pancreas_promoters.bed"
HUMAN_PANCREAS_ENHANCERS_BED="$CLASSIFIED_DIR/human/human_pancreas_enhancers.bed"
MOUSE_LIVER_PROMOTERS_BED="$CLASSIFIED_DIR/mouse/mouse_liver_promoters.bed"
MOUSE_LIVER_ENHANCERS_BED="$CLASSIFIED_DIR/mouse/mouse_liver_enhancers.bed"
MOUSE_PANCREAS_PROMOTERS_BED="$CLASSIFIED_DIR/mouse/mouse_pancreas_promoters.bed"
MOUSE_PANCREAS_ENHANCERS_BED="$CLASSIFIED_DIR/mouse/mouse_pancreas_enhancers.bed"
classify_peaks() { local input_peak_gz=$1; local promoter_def_bed=$2; local output_promoter_peaks_bed=$3; local output_enhancer_peaks_bed=$4; if [[ -s "$output_promoter_peaks_bed" && -s "$output_enhancer_peaks_bed" ]]; then echo "    Skipping classification for $(basename "$input_peak_gz"): Output files exist."; return; fi; local tmp_peak_bed="/tmp/$(basename "$input_peak_gz" .narrowPeak.gz)_$$_peaks.bed"; echo "    Classifying $(basename "$input_peak_gz") ..."; zcat "$input_peak_gz" | cut -f 1-3 > "$tmp_peak_bed"; bedtools intersect -a "$tmp_peak_bed" -b "$promoter_def_bed" -u > "$output_promoter_peaks_bed"; bedtools intersect -a "$tmp_peak_bed" -b "$promoter_def_bed" -v > "$output_enhancer_peaks_bed"; rm "$tmp_peak_bed"; echo "    Finished classifying $(basename "$input_peak_gz")"; }
classify_peaks "$HUMAN_LIVER_PEAKS" "$HUMAN_PROMOTERS_DEF_BED" "$HUMAN_LIVER_PROMOTERS_BED" "$HUMAN_LIVER_ENHANCERS_BED"
classify_peaks "$HUMAN_PANCREAS_PEAKS" "$HUMAN_PROMOTERS_DEF_BED" "$HUMAN_PANCREAS_PROMOTERS_BED" "$HUMAN_PANCREAS_ENHANCERS_BED"
classify_peaks "$MOUSE_LIVER_PEAKS" "$MOUSE_PROMOTERS_DEF_BED" "$MOUSE_LIVER_PROMOTERS_BED" "$MOUSE_LIVER_ENHANCERS_BED"
classify_peaks "$MOUSE_PANCREAS_PEAKS" "$MOUSE_PROMOTERS_DEF_BED" "$MOUSE_PANCREAS_PROMOTERS_BED" "$MOUSE_PANCREAS_ENHANCERS_BED"

# --- PART 4.3: Tissue Comparison of Enhancers ---
echo "  Performing tissue comparison on enhancer sets..."
HUMAN_ENHANCERS_TISSUE_SHARED="$COMP_DIR/tissue_comparison/human_enhancers_tissue_shared.bed"
HUMAN_LIVER_ENHANCERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/human_liver_enhancers_tissue_specific.bed"
HUMAN_PANCREAS_ENHANCERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/human_pancreas_enhancers_tissue_specific.bed"
MOUSE_ENHANCERS_TISSUE_SHARED="$COMP_DIR/tissue_comparison/mouse_enhancers_tissue_shared.bed"
MOUSE_LIVER_ENHANCERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/mouse_liver_enhancers_tissue_specific.bed"
MOUSE_PANCREAS_ENHANCERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/mouse_pancreas_enhancers_tissue_specific.bed"
bedtools intersect -a "$HUMAN_LIVER_ENHANCERS_BED" -b "$HUMAN_PANCREAS_ENHANCERS_BED" -u > "$HUMAN_ENHANCERS_TISSUE_SHARED"
bedtools intersect -a "$HUMAN_LIVER_ENHANCERS_BED" -b "$HUMAN_PANCREAS_ENHANCERS_BED" -v > "$HUMAN_LIVER_ENHANCERS_TISSUE_SPECIFIC"
bedtools intersect -a "$HUMAN_PANCREAS_ENHANCERS_BED" -b "$HUMAN_LIVER_ENHANCERS_BED" -v > "$HUMAN_PANCREAS_ENHANCERS_TISSUE_SPECIFIC"
bedtools intersect -a "$MOUSE_LIVER_ENHANCERS_BED" -b "$MOUSE_PANCREAS_ENHANCERS_BED" -u > "$MOUSE_ENHANCERS_TISSUE_SHARED"
bedtools intersect -a "$MOUSE_LIVER_ENHANCERS_BED" -b "$MOUSE_PANCREAS_ENHANCERS_BED" -v > "$MOUSE_LIVER_ENHANCERS_TISSUE_SPECIFIC"
bedtools intersect -a "$MOUSE_PANCREAS_ENHANCERS_BED" -b "$MOUSE_LIVER_ENHANCERS_BED" -v > "$MOUSE_PANCREAS_ENHANCERS_TISSUE_SPECIFIC"

# --- PART 4.4: Tissue Comparison of Promoters ---
echo "  Performing tissue comparison on promoter sets..."
HUMAN_PROMOTERS_TISSUE_SHARED="$COMP_DIR/tissue_comparison/human_promoters_tissue_shared.bed"
HUMAN_LIVER_PROMOTERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/human_liver_promoters_tissue_specific.bed"
HUMAN_PANCREAS_PROMOTERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/human_pancreas_promoters_tissue_specific.bed"
MOUSE_PROMOTERS_TISSUE_SHARED="$COMP_DIR/tissue_comparison/mouse_promoters_tissue_shared.bed"
MOUSE_LIVER_PROMOTERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/mouse_liver_promoters_tissue_specific.bed"
MOUSE_PANCREAS_PROMOTERS_TISSUE_SPECIFIC="$COMP_DIR/tissue_comparison/mouse_pancreas_promoters_tissue_specific.bed"
bedtools intersect -a "$HUMAN_LIVER_PROMOTERS_BED" -b "$HUMAN_PANCREAS_PROMOTERS_BED" -u > "$HUMAN_PROMOTERS_TISSUE_SHARED"
bedtools intersect -a "$HUMAN_LIVER_PROMOTERS_BED" -b "$HUMAN_PANCREAS_PROMOTERS_BED" -v > "$HUMAN_LIVER_PROMOTERS_TISSUE_SPECIFIC"
bedtools intersect -a "$HUMAN_PANCREAS_PROMOTERS_BED" -b "$HUMAN_LIVER_PROMOTERS_BED" -v > "$HUMAN_PANCREAS_PROMOTERS_TISSUE_SPECIFIC"
bedtools intersect -a "$MOUSE_LIVER_PROMOTERS_BED" -b "$MOUSE_PANCREAS_PROMOTERS_BED" -u > "$MOUSE_PROMOTERS_TISSUE_SHARED"
bedtools intersect -a "$MOUSE_LIVER_PROMOTERS_BED" -b "$MOUSE_PANCREAS_PROMOTERS_BED" -v > "$MOUSE_LIVER_PROMOTERS_TISSUE_SPECIFIC"
bedtools intersect -a "$MOUSE_PANCREAS_PROMOTERS_BED" -b "$MOUSE_LIVER_PROMOTERS_BED" -v > "$MOUSE_PANCREAS_PROMOTERS_TISSUE_SPECIFIC"

# --- PART 4.5: Species Comparison (Enhancers vs Mapped Peaks) ---
echo "  Performing species comparison (enhancers vs mapped peaks)..."
SKIP_SPECIES_STATS=0
if [[ ! -s "$MAPPED_ML_TO_HG_PEAKS" || ! -s "$MAPPED_MP_TO_HG_PEAKS" || ! -s "$MAPPED_HL_TO_MM_PEAKS" || ! -s "$MAPPED_HP_TO_MM_PEAKS" ]]; then echo "  Warning: One or more mapped peak BED files in $MAPPED_DIR are missing or empty. Skipping species comparison steps." >&2; SKIP_SPECIES_STATS=1; fi
if [[ "$SKIP_SPECIES_STATS" -eq 0 ]]; then
  HUMAN_LIVER_ENH_SHARED_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_liver_enhancers_shared_with_mapped_mouse_peaks_human_coords.bed"
  MOUSE_LIVER_ENH_SHARED_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_liver_enhancers_shared_with_mapped_human_peaks_mouse_coords.bed"
  HUMAN_LIVER_ENH_SPECIFIC_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_liver_enhancers_specific_vs_mapped_mouse_peaks_human_coords.bed"
  MOUSE_LIVER_ENH_SPECIFIC_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_liver_enhancers_specific_vs_mapped_human_peaks_mouse_coords.bed"
  HUMAN_PANCREAS_ENH_SHARED_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_pancreas_enhancers_shared_with_mapped_mouse_peaks_human_coords.bed"
  MOUSE_PANCREAS_ENH_SHARED_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_pancreas_enhancers_shared_with_mapped_human_peaks_mouse_coords.bed"
  HUMAN_PANCREAS_ENH_SPECIFIC_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_pancreas_enhancers_specific_vs_mapped_mouse_peaks_human_coords.bed"
  MOUSE_PANCREAS_ENH_SPECIFIC_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_pancreas_enhancers_specific_vs_mapped_human_peaks_mouse_coords.bed"
  bedtools intersect -a "$HUMAN_LIVER_ENHANCERS_BED" -b "$MAPPED_ML_TO_HG_PEAKS" -u > "$HUMAN_LIVER_ENH_SHARED_MAPPED_MOUSE_HG_COORDS"
  bedtools intersect -a "$MOUSE_LIVER_ENHANCERS_BED" -b "$MAPPED_HL_TO_MM_PEAKS" -u > "$MOUSE_LIVER_ENH_SHARED_MAPPED_HUMAN_MM_COORDS"
  bedtools intersect -a "$HUMAN_LIVER_ENHANCERS_BED" -b "$MAPPED_ML_TO_HG_PEAKS" -v > "$HUMAN_LIVER_ENH_SPECIFIC_MAPPED_MOUSE_HG_COORDS"
  bedtools intersect -a "$MOUSE_LIVER_ENHANCERS_BED" -b "$MAPPED_HL_TO_MM_PEAKS" -v > "$MOUSE_LIVER_ENH_SPECIFIC_MAPPED_HUMAN_MM_COORDS"
  bedtools intersect -a "$HUMAN_PANCREAS_ENHANCERS_BED" -b "$MAPPED_MP_TO_HG_PEAKS" -u > "$HUMAN_PANCREAS_ENH_SHARED_MAPPED_MOUSE_HG_COORDS"
  bedtools intersect -a "$MOUSE_PANCREAS_ENHANCERS_BED" -b "$MAPPED_HP_TO_MM_PEAKS" -u > "$MOUSE_PANCREAS_ENH_SHARED_MAPPED_HUMAN_MM_COORDS"
  bedtools intersect -a "$HUMAN_PANCREAS_ENHANCERS_BED" -b "$MAPPED_MP_TO_HG_PEAKS" -v > "$HUMAN_PANCREAS_ENH_SPECIFIC_MAPPED_MOUSE_HG_COORDS"
  bedtools intersect -a "$MOUSE_PANCREAS_ENHANCERS_BED" -b "$MAPPED_HP_TO_MM_PEAKS" -v > "$MOUSE_PANCREAS_ENH_SPECIFIC_MAPPED_HUMAN_MM_COORDS"
else echo "  Skipping species comparison for enhancers."; fi

# --- PART 4.6: Species Comparison (Promoters vs Mapped Peaks) ---
echo "  Performing species comparison (promoters vs mapped peaks)..."
if [[ "$SKIP_SPECIES_STATS" -eq 0 ]]; then
  HUMAN_LIVER_PROM_SHARED_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_liver_promoters_shared_with_mapped_mouse_peaks_human_coords.bed"
  MOUSE_LIVER_PROM_SHARED_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_liver_promoters_shared_with_mapped_human_peaks_mouse_coords.bed"
  HUMAN_LIVER_PROM_SPECIFIC_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_liver_promoters_specific_vs_mapped_mouse_peaks_human_coords.bed"
  MOUSE_LIVER_PROM_SPECIFIC_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_liver_promoters_specific_vs_mapped_human_peaks_mouse_coords.bed"
  HUMAN_PANCREAS_PROM_SHARED_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_pancreas_promoters_shared_with_mapped_mouse_peaks_human_coords.bed"
  MOUSE_PANCREAS_PROM_SHARED_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_pancreas_promoters_shared_with_mapped_human_peaks_mouse_coords.bed"
  HUMAN_PANCREAS_PROM_SPECIFIC_MAPPED_MOUSE_HG_COORDS="$COMP_DIR/species_comparison/human_pancreas_promoters_specific_vs_mapped_mouse_peaks_human_coords.bed"
  MOUSE_PANCREAS_PROM_SPECIFIC_MAPPED_HUMAN_MM_COORDS="$COMP_DIR/species_comparison/mouse_pancreas_promoters_specific_vs_mapped_human_peaks_mouse_coords.bed"
  bedtools intersect -a "$HUMAN_LIVER_PROMOTERS_BED" -b "$MAPPED_ML_TO_HG_PEAKS" -u > "$HUMAN_LIVER_PROM_SHARED_MAPPED_MOUSE_HG_COORDS"
  bedtools intersect -a "$MOUSE_LIVER_PROMOTERS_BED" -b "$MAPPED_HL_TO_MM_PEAKS" -u > "$MOUSE_LIVER_PROM_SHARED_MAPPED_HUMAN_MM_COORDS"
  bedtools intersect -a "$HUMAN_LIVER_PROMOTERS_BED" -b "$MAPPED_ML_TO_HG_PEAKS" -v > "$HUMAN_LIVER_PROM_SPECIFIC_MAPPED_MOUSE_HG_COORDS"
  bedtools intersect -a "$MOUSE_LIVER_PROMOTERS_BED" -b "$MAPPED_HL_TO_MM_PEAKS" -v > "$MOUSE_LIVER_PROM_SPECIFIC_MAPPED_HUMAN_MM_COORDS"
  bedtools intersect -a "$HUMAN_PANCREAS_PROMOTERS_BED" -b "$MAPPED_MP_TO_HG_PEAKS" -u > "$HUMAN_PANCREAS_PROM_SHARED_MAPPED_MOUSE_HG_COORDS"
  bedtools intersect -a "$MOUSE_PANCREAS_PROMOTERS_BED" -b "$MAPPED_HP_TO_MM_PEAKS" -u > "$MOUSE_PANCREAS_PROM_SHARED_MAPPED_HUMAN_MM_COORDS"
  bedtools intersect -a "$HUMAN_PANCREAS_PROMOTERS_BED" -b "$MAPPED_MP_TO_HG_PEAKS" -v > "$HUMAN_PANCREAS_PROM_SPECIFIC_MAPPED_MOUSE_HG_COORDS"
  bedtools intersect -a "$MOUSE_PANCREAS_PROMOTERS_BED" -b "$MAPPED_HP_TO_MM_PEAKS" -v > "$MOUSE_PANCREAS_PROM_SPECIFIC_MAPPED_HUMAN_MM_COORDS"
else echo "  Skipping species comparison for promoters."; fi

# --- PART 4.7: Calculate Statistics & Percentages ---
echo "  Calculating statistics and percentages..."
STATS_FILE="$ANALYSIS_OUTPUT_DIR/classification_and_comparison_summary.txt"
calc_pct() { local num_str=$(echo "$1" | awk '{$1=$1;print}'); local denom_str=$(echo "$2" | awk '{$1=$1;print}'); local num=${num_str:-0}; local denom=${denom_str:-0}; if [[ "$denom" -eq 0 ]]; then echo "0.00"; else echo "scale=2; ($num / $denom) * 100" | bc; fi; }
count_lines() { cat "$1" 2>/dev/null | wc -l; }
HLP_C=$(count_lines "$HUMAN_LIVER_PROMOTERS_BED"); HLE_C=$(count_lines "$HUMAN_LIVER_ENHANCERS_BED")
HPP_C=$(count_lines "$HUMAN_PANCREAS_PROMOTERS_BED"); HPE_C=$(count_lines "$HUMAN_PANCREAS_ENHANCERS_BED")
MLP_C=$(count_lines "$MOUSE_LIVER_PROMOTERS_BED"); MLE_C=$(count_lines "$MOUSE_LIVER_ENHANCERS_BED")
MPP_C=$(count_lines "$MOUSE_PANCREAS_PROMOTERS_BED"); MPE_C=$(count_lines "$MOUSE_PANCREAS_ENHANCERS_BED")
HETS_C=$(count_lines "$HUMAN_ENHANCERS_TISSUE_SHARED")
HPTS_C=$(count_lines "$HUMAN_PROMOTERS_TISSUE_SHARED")
METS_C=$(count_lines "$MOUSE_ENHANCERS_TISSUE_SHARED")
MPTS_C=$(count_lines "$MOUSE_PROMOTERS_TISSUE_SHARED")

HLES_C=0; MLES_C=0; HPES_C=0; MPES_C=0; HPLS_C=0; MPLS_C=0; HPPS_C=0; MPPS_C=0
if [[ "$SKIP_SPECIES_STATS" -eq 0 ]]; then
  HLES_C=$(count_lines "$HUMAN_LIVER_ENH_SHARED_MAPPED_MOUSE_HG_COORDS")
  MLES_C=$(count_lines "$MOUSE_LIVER_ENH_SHARED_MAPPED_HUMAN_MM_COORDS")
  HPES_C=$(count_lines "$HUMAN_PANCREAS_ENH_SHARED_MAPPED_MOUSE_HG_COORDS")
  MPES_C=$(count_lines "$MOUSE_PANCREAS_ENH_SHARED_MAPPED_HUMAN_MM_COORDS")
  HPLS_C=$(count_lines "$HUMAN_LIVER_PROM_SHARED_MAPPED_MOUSE_HG_COORDS")
  MPLS_C=$(count_lines "$MOUSE_LIVER_PROM_SHARED_MAPPED_HUMAN_MM_COORDS")
  HPPS_C=$(count_lines "$HUMAN_PANCREAS_PROM_SHARED_MAPPED_MOUSE_HG_COORDS")
  MPPS_C=$(count_lines "$MOUSE_PANCREAS_PROM_SHARED_MAPPED_HUMAN_MM_COORDS")
fi

HL_ENH_SHARED_PCT=$(calc_pct $HETS_C $HLE_C); HP_ENH_SHARED_PCT=$(calc_pct $HETS_C $HPE_C)
HL_PROM_SHARED_PCT=$(calc_pct $HPTS_C $HLP_C); HP_PROM_SHARED_PCT=$(calc_pct $HPTS_C $HPP_C)
ML_ENH_SHARED_PCT=$(calc_pct $METS_C $MLE_C); MP_ENH_SHARED_PCT=$(calc_pct $METS_C $MPE_C)
ML_PROM_SHARED_PCT=$(calc_pct $MPTS_C $MLP_C); MP_PROM_SHARED_PCT=$(calc_pct $MPTS_C $MPP_C)
HL_ENH_SHARED_SPECIES_PCT=$(calc_pct $HLES_C $HLE_C); HP_ENH_SHARED_SPECIES_PCT=$(calc_pct $HPES_C $HPE_C)
HL_PROM_SHARED_SPECIES_PCT=$(calc_pct $HPLS_C $HLP_C); HP_PROM_SHARED_SPECIES_PCT=$(calc_pct $HPPS_C $HPP_C)
ML_ENH_SHARED_SPECIES_PCT=$(calc_pct $MLES_C $MLE_C); MP_ENH_SHARED_SPECIES_PCT=$(calc_pct $MPES_C $MPE_C)
ML_PROM_SHARED_SPECIES_PCT=$(calc_pct $MPLS_C $MLP_C); MP_PROM_SHARED_SPECIES_PCT=$(calc_pct $MPPS_C $MPP_C)

echo "===== Classification and Comparison Summary =====" > "$STATS_FILE"
echo "Timestamp: $(date)" >> "$STATS_FILE"
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

echo "  Classification and comparison statistics written to $STATS_FILE"

echo "===== Completed Step 4 ====="