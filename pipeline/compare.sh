#!/bin/bash
# compare.sh - Performs initial cross-species and tissue comparison using bedtools.
#              Processes HALPER narrowPeak output and compares with original peaks.
#
# Usage: ./compare.sh <output_dir> <human_liver_peaks_gz> <human_pancreas_peaks_gz> <mouse_liver_peaks_gz> <mouse_pancreas_peaks_gz> <human_liver_to_mouse_halper_np> <human_pancreas_to_mouse_halper_np> <mouse_liver_to_human_halper_np> <mouse_pancreas_to_human_halper_np>
#
# Inputs:
#   Arg 1:  <output_dir>                        - Base directory for pipeline results.
#   Arg 2:  <human_liver_peaks_gz>            - Gzipped narrowPeak: Original Human Liver peaks.
#   Arg 3:  <human_pancreas_peaks_gz>         - Gzipped narrowPeak: Original Human Pancreas peaks.
#   Arg 4:  <mouse_liver_peaks_gz>            - Gzipped narrowPeak: Original Mouse Liver peaks.
#   Arg 5:  <mouse_pancreas_peaks_gz>         - Gzipped narrowPeak: Original Mouse Pancreas peaks.
#   Arg 6:  <human_liver_to_mouse_halper_np>  - Gzipped narrowPeak: HALPER output (Human Liver -> Mouse).
#   Arg 7:  <human_pancreas_to_mouse_halper_np> - Gzipped narrowPeak: HALPER output (Human Pancreas -> Mouse).
#   Arg 8:  <mouse_liver_to_human_halper_np>  - Gzipped narrowPeak: HALPER output (Mouse Liver -> Human).
#   Arg 9:  <mouse_pancreas_to_human_halper_np> - Gzipped narrowPeak: HALPER output (Mouse Pancreas -> Human).
#
# Outputs:
#   Directory: <output_dir>/initial_comparisons/tissue_comparison/          - Contains BED files from tissue comparisons.
#   File:      .../human_shared_between_liver_and_pancreas.bed - Human peaks overlapping between tissues.
#   File:      .../human_liver_specific.bed                    - Human peaks unique to liver.
#   File:      .../human_pancreas_specific.bed                 - Human peaks unique to pancreas.
#   File:      .../mouse_shared_between_liver_and_pancreas.bed - Mouse peaks overlapping between tissues.
#   File:      .../mouse_liver_specific.bed                    - Mouse peaks unique to liver.
#   File:      .../mouse_pancreas_specific.bed                 - Mouse peaks unique to pancreas.
#   File:      <output_dir>/initial_comparisons/step2_tissue_comparison_summary.txt - Summary statistics for tissue comparisons.
#
#   Directory: <output_dir>/initial_comparisons/cross_species/            - Contains BED files from cross-species comparisons.
#   Directory: .../processed_files/                                       - Contains BED files derived from HALPER narrowPeak inputs.
#   File:      .../human_liver_conserved_in_mouse_liver.bed         - Mapped Human liver peaks overlapping Mouse liver peaks.
#   File:      .../human_liver_conserved_in_mouse_pancreas.bed      - Mapped Human liver peaks overlapping Mouse pancreas peaks.
#   File:      .../human_pancreas_conserved_in_mouse_pancreas.bed   - Mapped Human pancreas peaks overlapping Mouse pancreas peaks.
#   File:      .../human_pancreas_conserved_in_mouse_liver.bed      - Mapped Human pancreas peaks overlapping Mouse liver peaks.
#   File:      .../mouse_liver_conserved_in_human_liver.bed         - Mapped Mouse liver peaks overlapping Human liver peaks.
#   File:      .../mouse_liver_conserved_in_human_pancreas.bed      - Mapped Mouse liver peaks overlapping Human pancreas peaks.
#   File:      .../mouse_pancreas_conserved_in_human_pancreas.bed   - Mapped Mouse pancreas peaks overlapping Human pancreas peaks.
#   File:      .../mouse_pancreas_conserved_in_human_liver.bed      - Mapped Mouse pancreas peaks overlapping Human liver peaks.
#   File:      .../human_liver_not_conserved_in_mouse.bed           - Mapped Human liver peaks not overlapping any Mouse peaks.
#   File:      .../human_pancreas_not_conserved_in_mouse.bed        - Mapped Human pancreas peaks not overlapping any Mouse peaks.
#   File:      .../mouse_liver_not_conserved_in_human.bed           - Mapped Mouse liver peaks not overlapping any Human peaks.
#   File:      .../mouse_pancreas_not_conserved_in_human.bed        - Mapped Mouse pancreas peaks not overlapping any Human peaks.
#   File:      <output_dir>/initial_comparisons/step2_cross_species_conservation_summary.txt - Summary statistics for cross-species comparisons.

set -e

# Parse command-line arguments
if [ $# -lt 9 ]; then
    echo "Usage: $0 <output_dir> <human_liver_peaks_gz> <human_pancreas_peaks_gz> <mouse_liver_peaks_gz> <mouse_pancreas_peaks_gz> <human_liver_to_mouse_halper_np> <human_pancreas_to_mouse_halper_np> <mouse_liver_to_human_halper_np> <mouse_pancreas_to_human_halper_np>"
    exit 1
fi

OUTPUT_DIR=$1
HUMAN_LIVER_PEAKS=$2
HUMAN_PANCREAS_PEAKS=$3
MOUSE_LIVER_PEAKS=$4
MOUSE_PANCREAS_PEAKS=$5
HUMAN_LIVER_TO_MOUSE_HALPER=$6
HUMAN_PANCREAS_TO_MOUSE_HALPER=$7
MOUSE_LIVER_TO_HUMAN_HALPER=$8
MOUSE_PANCREAS_TO_HUMAN_HALPER=$9

# Define output subdirectories
STEP2_OUTPUT_DIR="$OUTPUT_DIR/initial_comparisons"
CROSS_SPECIES_DIR="$STEP2_OUTPUT_DIR/cross_species"
TISSUE_COMP_DIR="$STEP2_OUTPUT_DIR/tissue_comparison"
PROCESSED_DIR="$CROSS_SPECIES_DIR/processed_files"
mkdir -p "$STEP2_OUTPUT_DIR" "$CROSS_SPECIES_DIR" "$TISSUE_COMP_DIR" "$PROCESSED_DIR"

# Log input files
echo "  Inputs for Step 2:"
echo "    Human liver peaks: $HUMAN_LIVER_PEAKS"
echo "    Human pancreas peaks: $HUMAN_PANCREAS_PEAKS"
echo "    Mouse liver peaks: $MOUSE_LIVER_PEAKS"
echo "    Mouse pancreas peaks: $MOUSE_PANCREAS_PEAKS"
echo "    Human liver to mouse HALPER NP: $HUMAN_LIVER_TO_MOUSE_HALPER"
echo "    Human pancreas to mouse HALPER NP: $HUMAN_PANCREAS_TO_MOUSE_HALPER"
echo "    Mouse liver to human HALPER NP: $MOUSE_LIVER_TO_HUMAN_HALPER"
echo "    Mouse pancreas to human HALPER NP: $MOUSE_PANCREAS_TO_HUMAN_HALPER"

# Check for required input files
all_step2_inputs_found=true
for file in "$HUMAN_LIVER_PEAKS" "$HUMAN_PANCREAS_PEAKS" "$MOUSE_LIVER_PEAKS" "$MOUSE_PANCREAS_PEAKS" \
            "$HUMAN_LIVER_TO_MOUSE_HALPER" "$HUMAN_PANCREAS_TO_MOUSE_HALPER" \
            "$MOUSE_LIVER_TO_HUMAN_HALPER" "$MOUSE_PANCREAS_TO_HUMAN_HALPER"; do
    if [ ! -f "$file" ]; then
        echo "Error: Input file $file required for Step 2 does not exist" >&2
        all_step2_inputs_found=false
    fi
done
if ! $all_step2_inputs_found; then exit 1; fi

# Function to safely calculate percentages
calc_pct_step2() {
    local num_str=$(echo "$1" | awk '{$1=$1;print}')
    local denom_str=$(echo "$2" | awk '{$1=$1;print}')
    local num=${num_str:-0}
    local denom=${denom_str:-0}
    if [[ "$denom" -eq 0 ]]; then echo "0.00"; else echo "scale=2; ($num / $denom) * 100" | bc; fi
}

echo "===== STEP 2: Running Initial Bedtools Comparisons ====="

# --- PART 2.1: WITHIN-SPECIES TISSUE COMPARISON (using original peaks) ---
echo "  Starting Tissue Comparison Analysis (using original peaks)..."

# Convert narrowPeak files to BED format
echo "  Converting narrowPeak files to BED format..."
HUMAN_LIVER_BED="$TISSUE_COMP_DIR/human_liver.bed"
HUMAN_PANCREAS_BED="$TISSUE_COMP_DIR/human_pancreas.bed"
MOUSE_LIVER_BED="$TISSUE_COMP_DIR/mouse_liver.bed"
MOUSE_PANCREAS_BED="$TISSUE_COMP_DIR/mouse_pancreas.bed"
zcat "$HUMAN_LIVER_PEAKS" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,".",$7}' > "$HUMAN_LIVER_BED" # Keep score column ($7)
zcat "$HUMAN_PANCREAS_PEAKS" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,".",$7}' > "$HUMAN_PANCREAS_BED"
zcat "$MOUSE_LIVER_PEAKS" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,".",$7}' > "$MOUSE_LIVER_BED"
zcat "$MOUSE_PANCREAS_PEAKS" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,".",$7}' > "$MOUSE_PANCREAS_BED"

# Human tissue comparison
echo "  Finding shared and tissue-specific regions in human..."
HUMAN_SHARED_STEP2="$TISSUE_COMP_DIR/human_shared_between_liver_and_pancreas.bed"
HUMAN_LIVER_SPECIFIC_STEP2="$TISSUE_COMP_DIR/human_liver_specific.bed"
HUMAN_PANCREAS_SPECIFIC_STEP2="$TISSUE_COMP_DIR/human_pancreas_specific.bed"
bedtools intersect -a "$HUMAN_LIVER_BED" -b "$HUMAN_PANCREAS_BED" -u > "$HUMAN_SHARED_STEP2"
bedtools intersect -a "$HUMAN_LIVER_BED" -b "$HUMAN_PANCREAS_BED" -v > "$HUMAN_LIVER_SPECIFIC_STEP2"
bedtools intersect -a "$HUMAN_PANCREAS_BED" -b "$HUMAN_LIVER_BED" -v > "$HUMAN_PANCREAS_SPECIFIC_STEP2"

# Mouse tissue comparison
echo "  Finding shared and tissue-specific regions in mouse..."
MOUSE_SHARED_STEP2="$TISSUE_COMP_DIR/mouse_shared_between_liver_and_pancreas.bed"
MOUSE_LIVER_SPECIFIC_STEP2="$TISSUE_COMP_DIR/mouse_liver_specific.bed"
MOUSE_PANCREAS_SPECIFIC_STEP2="$TISSUE_COMP_DIR/mouse_pancreas_specific.bed"
bedtools intersect -a "$MOUSE_LIVER_BED" -b "$MOUSE_PANCREAS_BED" -u > "$MOUSE_SHARED_STEP2"
bedtools intersect -a "$MOUSE_LIVER_BED" -b "$MOUSE_PANCREAS_BED" -v > "$MOUSE_LIVER_SPECIFIC_STEP2"
bedtools intersect -a "$MOUSE_PANCREAS_BED" -b "$MOUSE_LIVER_BED" -v > "$MOUSE_PANCREAS_SPECIFIC_STEP2"

# Calculate tissue sharing statistics
echo "  Calculating tissue sharing statistics..."
HL_COUNT_S2=$(wc -l < "$HUMAN_LIVER_BED"); HP_COUNT_S2=$(wc -l < "$HUMAN_PANCREAS_BED")
HS_COUNT_S2=$(wc -l < "$HUMAN_SHARED_STEP2"); HLS_COUNT_S2=$(wc -l < "$HUMAN_LIVER_SPECIFIC_STEP2"); HPS_COUNT_S2=$(wc -l < "$HUMAN_PANCREAS_SPECIFIC_STEP2")
ML_COUNT_S2=$(wc -l < "$MOUSE_LIVER_BED"); MP_COUNT_S2=$(wc -l < "$MOUSE_PANCREAS_BED")
MS_COUNT_S2=$(wc -l < "$MOUSE_SHARED_STEP2"); MLS_COUNT_S2=$(wc -l < "$MOUSE_LIVER_SPECIFIC_STEP2"); MPS_COUNT_S2=$(wc -l < "$MOUSE_PANCREAS_SPECIFIC_STEP2")
# Calculate percentages
HL_SHARED_PCT_S2=$(calc_pct_step2 $HS_COUNT_S2 $HL_COUNT_S2); HP_SHARED_PCT_S2=$(calc_pct_step2 $HS_COUNT_S2 $HP_COUNT_S2)
HLS_PCT_S2=$(calc_pct_step2 $HLS_COUNT_S2 $HL_COUNT_S2); HPS_PCT_S2=$(calc_pct_step2 $HPS_COUNT_S2 $HP_COUNT_S2)
ML_SHARED_PCT_S2=$(calc_pct_step2 $MS_COUNT_S2 $ML_COUNT_S2); MP_SHARED_PCT_S2=$(calc_pct_step2 $MS_COUNT_S2 $MP_COUNT_S2)
MLS_PCT_S2=$(calc_pct_step2 $MLS_COUNT_S2 $ML_COUNT_S2); MPS_PCT_S2=$(calc_pct_step2 $MPS_COUNT_S2 $MP_COUNT_S2)

# Create tissue comparison summary
TISSUE_SUMMARY_FILE="$STEP2_OUTPUT_DIR/step2_tissue_comparison_summary.txt"
echo "===== STEP 2: TISSUE COMPARISON SUMMARY =====" > "$TISSUE_SUMMARY_FILE"; echo "" >> "$TISSUE_SUMMARY_FILE"
echo "HUMAN TISSUE COMPARISON" >> "$TISSUE_SUMMARY_FILE"; echo "-----------------------------" >> "$TISSUE_SUMMARY_FILE"
echo "Total liver regions: $HL_COUNT_S2" >> "$TISSUE_SUMMARY_FILE"; echo "Total pancreas regions: $HP_COUNT_S2" >> "$TISSUE_SUMMARY_FILE"
echo "Regions shared between tissues: $HS_COUNT_S2" >> "$TISSUE_SUMMARY_FILE"
echo "Liver-specific regions: $HLS_COUNT_S2 - ${HLS_PCT_S2}%" >> "$TISSUE_SUMMARY_FILE"; echo "Pancreas-specific regions: $HPS_COUNT_S2 - ${HPS_PCT_S2}%" >> "$TISSUE_SUMMARY_FILE"
echo "Liver sharing with pancreas: ${HL_SHARED_PCT_S2}%" >> "$TISSUE_SUMMARY_FILE"; echo "Pancreas sharing with liver: ${HP_SHARED_PCT_S2}%" >> "$TISSUE_SUMMARY_FILE"; echo "" >> "$TISSUE_SUMMARY_FILE"
echo "MOUSE TISSUE COMPARISON" >> "$TISSUE_SUMMARY_FILE"; echo "-----------------------------" >> "$TISSUE_SUMMARY_FILE"
echo "Total liver regions: $ML_COUNT_S2" >> "$TISSUE_SUMMARY_FILE"; echo "Total pancreas regions: $MP_COUNT_S2" >> "$TISSUE_SUMMARY_FILE"
echo "Regions shared between tissues: $MS_COUNT_S2" >> "$TISSUE_SUMMARY_FILE"
echo "Liver-specific regions: $MLS_COUNT_S2 - ${MLS_PCT_S2}%" >> "$TISSUE_SUMMARY_FILE"; echo "Pancreas-specific regions: $MPS_COUNT_S2 - ${MPS_PCT_S2}%" >> "$TISSUE_SUMMARY_FILE"
echo "Liver sharing with pancreas: ${ML_SHARED_PCT_S2}%" >> "$TISSUE_SUMMARY_FILE"; echo "Pancreas sharing with liver: ${MP_SHARED_PCT_S2}%" >> "$TISSUE_SUMMARY_FILE"; echo "" >> "$TISSUE_SUMMARY_FILE"
echo "  Tissue comparison complete. Results in $TISSUE_COMP_DIR. Summary in $TISSUE_SUMMARY_FILE"

# --- PART 2.2: CROSS-SPECIES CONSERVATION ANALYSIS (using original peaks and processed HALPER output) ---
echo "  Starting Cross-Species Conservation Analysis..."

# Extract properly formatted coordinates from HALPER narrowPeak output
echo "  Processing HALPER output files to BED format..."
PROCESSED_HL_TO_MM="$PROCESSED_DIR/human_liver_to_mouse.bed"
PROCESSED_HP_TO_MM="$PROCESSED_DIR/human_pancreas_to_mouse.bed"
PROCESSED_ML_TO_HG="$PROCESSED_DIR/mouse_liver_to_human.bed"
PROCESSED_MP_TO_HG="$PROCESSED_DIR/mouse_pancreas_to_human.bed"

# Function to process HALPER narrowPeak files and extract mapped coordinates to BED
process_halper_file() {
    local input_np_gz=$1; local output_bed=$2; local line_count
    echo "    Processing $(basename "$input_np_gz") -> $(basename "$output_bed")"
    zcat "$input_np_gz" | awk 'BEGIN{OFS="\t"} {
        # HALPER narrowPeak Col 4 format: chr_q:start_q-end_q:summit_q
        # We want the mapped coordinates (target genome) which are in Cols 1-3
        # We store original query coords and summit in Col 4 for reference if needed
        # We extract the p-value (-log10) from Col 8 ($8) as the score for Col 5
        print $1, $2, $3, $4, $8, "." # Output: chr_t, start_t, end_t, name (query_coords), score (p-value), strand (.)
    }' > "$output_bed"
    line_count=$(wc -l < "$output_bed")
    echo "      Processed $input_np_gz: $line_count regions written to $output_bed"
}

process_halper_file "$HUMAN_LIVER_TO_MOUSE_HALPER" "$PROCESSED_HL_TO_MM"
process_halper_file "$HUMAN_PANCREAS_TO_MOUSE_HALPER" "$PROCESSED_HP_TO_MM"
process_halper_file "$MOUSE_LIVER_TO_HUMAN_HALPER" "$PROCESSED_ML_TO_HG"
process_halper_file "$MOUSE_PANCREAS_TO_HUMAN_HALPER" "$PROCESSED_MP_TO_HG"

# Find conserved regions using the processed HALPER BED output vs original peak BEDs
echo "  Finding conserved regions across species..."
# Human liver regions conserved in mouse liver/pancreas
HL_CONS_ML="$CROSS_SPECIES_DIR/human_liver_conserved_in_mouse_liver.bed" # Mapped coords (mouse) that overlap original mouse liver peaks
HL_CONS_MP="$CROSS_SPECIES_DIR/human_liver_conserved_in_mouse_pancreas.bed"
bedtools intersect -a "$PROCESSED_HL_TO_MM" -b "$MOUSE_LIVER_BED" -u > "$HL_CONS_ML"
bedtools intersect -a "$PROCESSED_HL_TO_MM" -b "$MOUSE_PANCREAS_BED" -u > "$HL_CONS_MP"
# Human pancreas regions conserved in mouse liver/pancreas
HP_CONS_MP="$CROSS_SPECIES_DIR/human_pancreas_conserved_in_mouse_pancreas.bed"
HP_CONS_ML="$CROSS_SPECIES_DIR/human_pancreas_conserved_in_mouse_liver.bed"
bedtools intersect -a "$PROCESSED_HP_TO_MM" -b "$MOUSE_PANCREAS_BED" -u > "$HP_CONS_MP"
bedtools intersect -a "$PROCESSED_HP_TO_MM" -b "$MOUSE_LIVER_BED" -u > "$HP_CONS_ML"
# Mouse liver regions conserved in human liver/pancreas
ML_CONS_HL="$CROSS_SPECIES_DIR/mouse_liver_conserved_in_human_liver.bed"
ML_CONS_HP="$CROSS_SPECIES_DIR/mouse_liver_conserved_in_human_pancreas.bed"
bedtools intersect -a "$PROCESSED_ML_TO_HG" -b "$HUMAN_LIVER_BED" -u > "$ML_CONS_HL"
bedtools intersect -a "$PROCESSED_ML_TO_HG" -b "$HUMAN_PANCREAS_BED" -u > "$ML_CONS_HP"
# Mouse pancreas regions conserved in human liver/pancreas
MP_CONS_HP="$CROSS_SPECIES_DIR/mouse_pancreas_conserved_in_human_pancreas.bed"
MP_CONS_HL="$CROSS_SPECIES_DIR/mouse_pancreas_conserved_in_human_liver.bed"
bedtools intersect -a "$PROCESSED_MP_TO_HG" -b "$HUMAN_PANCREAS_BED" -u > "$MP_CONS_HP"
bedtools intersect -a "$PROCESSED_MP_TO_HG" -b "$HUMAN_LIVER_BED" -u > "$MP_CONS_HL"

# Find non-conserved/unique regions
echo "  Finding species-specific unique regions (mapped regions not overlapping any peak in target)..."
# Create combined BED files for *all* original peaks in each target species
MOUSE_ALL_PEAKS_BED="$PROCESSED_DIR/mouse_all_tissues_merged.bed"
HUMAN_ALL_PEAKS_BED="$PROCESSED_DIR/human_all_tissues_merged.bed"
cat "$MOUSE_LIVER_BED" "$MOUSE_PANCREAS_BED" | sort -k1,1 -k2,2n | bedtools merge -i - > "$MOUSE_ALL_PEAKS_BED"
cat "$HUMAN_LIVER_BED" "$HUMAN_PANCREAS_BED" | sort -k1,1 -k2,2n | bedtools merge -i - > "$HUMAN_ALL_PEAKS_BED"
# Human mapped regions not conserved in mouse
HL_NOT_CONS_S2="$CROSS_SPECIES_DIR/human_liver_not_conserved_in_mouse.bed" # Mapped coords (mouse) that *don't* overlap any mouse peak
HP_NOT_CONS_S2="$CROSS_SPECIES_DIR/human_pancreas_not_conserved_in_mouse.bed"
bedtools intersect -a "$PROCESSED_HL_TO_MM" -b "$MOUSE_ALL_PEAKS_BED" -v > "$HL_NOT_CONS_S2"
bedtools intersect -a "$PROCESSED_HP_TO_MM" -b "$MOUSE_ALL_PEAKS_BED" -v > "$HP_NOT_CONS_S2"
# Mouse mapped regions not conserved in human
ML_NOT_CONS_S2="$CROSS_SPECIES_DIR/mouse_liver_not_conserved_in_human.bed"
MP_NOT_CONS_S2="$CROSS_SPECIES_DIR/mouse_pancreas_not_conserved_in_human.bed"
bedtools intersect -a "$PROCESSED_ML_TO_HG" -b "$HUMAN_ALL_PEAKS_BED" -v > "$ML_NOT_CONS_S2"
bedtools intersect -a "$PROCESSED_MP_TO_HG" -b "$HUMAN_ALL_PEAKS_BED" -v > "$MP_NOT_CONS_S2"

# Calculate cross-species conservation statistics
echo "  Calculating cross-species conservation statistics..."
# Counts of mapped regions
HL_MAP_C_S2=$(wc -l < "$PROCESSED_HL_TO_MM"); HP_MAP_C_S2=$(wc -l < "$PROCESSED_HP_TO_MM")
ML_MAP_C_S2=$(wc -l < "$PROCESSED_ML_TO_HG"); MP_MAP_C_S2=$(wc -l < "$PROCESSED_MP_TO_HG")
# Counts of conserved regions
HL_C_ML_C=$(wc -l < "$HL_CONS_ML"); HL_C_MP_C=$(wc -l < "$HL_CONS_MP")
HP_C_MP_C=$(wc -l < "$HP_CONS_MP"); HP_C_ML_C=$(wc -l < "$HP_CONS_ML")
ML_C_HL_C=$(wc -l < "$ML_CONS_HL"); ML_C_HP_C=$(wc -l < "$ML_CONS_HP")
MP_C_HP_C=$(wc -l < "$MP_CONS_HP"); MP_C_HL_C=$(wc -l < "$MP_CONS_HL")
# Counts of non-conserved regions
HL_NC_C_S2=$(wc -l < "$HL_NOT_CONS_S2"); HP_NC_C_S2=$(wc -l < "$HP_NOT_CONS_S2")
ML_NC_C_S2=$(wc -l < "$ML_NOT_CONS_S2"); MP_NC_C_S2=$(wc -l < "$MP_NOT_CONS_S2")
# Calculate percentages
HL_LIFT_PCT_S2=$(calc_pct_step2 $HL_MAP_C_S2 $HL_COUNT_S2); HP_LIFT_PCT_S2=$(calc_pct_step2 $HP_MAP_C_S2 $HP_COUNT_S2)
ML_LIFT_PCT_S2=$(calc_pct_step2 $ML_MAP_C_S2 $ML_COUNT_S2); MP_LIFT_PCT_S2=$(calc_pct_step2 $MP_MAP_C_S2 $MP_COUNT_S2)
HL_C_ML_PCT_S2=$(calc_pct_step2 $HL_C_ML_C $HL_MAP_C_S2); HL_C_MP_PCT_S2=$(calc_pct_step2 $HL_C_MP_C $HL_MAP_C_S2)
HP_C_MP_PCT_S2=$(calc_pct_step2 $HP_C_MP_C $HP_MAP_C_S2); HP_C_ML_PCT_S2=$(calc_pct_step2 $HP_C_ML_C $HP_MAP_C_S2)
ML_C_HL_PCT_S2=$(calc_pct_step2 $ML_C_HL_C $ML_MAP_C_S2); ML_C_HP_PCT_S2=$(calc_pct_step2 $ML_C_HP_C $ML_MAP_C_S2)
MP_C_HP_PCT_S2=$(calc_pct_step2 $MP_C_HP_C $MP_MAP_C_S2); MP_C_HL_PCT_S2=$(calc_pct_step2 $MP_C_HL_C $MP_MAP_C_S2)
HL_UNIQUE_PCT_S2=$(calc_pct_step2 $HL_NC_C_S2 $HL_MAP_C_S2); HP_UNIQUE_PCT_S2=$(calc_pct_step2 $HP_NC_C_S2 $HP_MAP_C_S2)
ML_UNIQUE_PCT_S2=$(calc_pct_step2 $ML_NC_C_S2 $ML_MAP_C_S2); MP_UNIQUE_PCT_S2=$(calc_pct_step2 $MP_NC_C_S2 $MP_MAP_C_S2)

# Create cross-species conservation summary
CROSS_SPECIES_SUMMARY_FILE="$STEP2_OUTPUT_DIR/step2_cross_species_conservation_summary.txt"
echo "===== STEP 2: CROSS-SPECIES CONSERVATION SUMMARY =====" > "$CROSS_SPECIES_SUMMARY_FILE"; echo "" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "HUMAN TO MOUSE CONSERVATION STATISTICS" >> "$CROSS_SPECIES_SUMMARY_FILE"; echo "-----------------------------" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Total human liver regions: $HL_COUNT_S2" >> "$CROSS_SPECIES_SUMMARY_FILE"; echo "Total human pancreas regions: $HP_COUNT_S2" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human liver regions mapped to mouse: $HL_MAP_C_S2 - ${HL_LIFT_PCT_S2}%" >> "$CROSS_SPECIES_SUMMARY_FILE"; echo "Human pancreas regions mapped to mouse: $HP_MAP_C_S2 - ${HP_LIFT_PCT_S2}%" >> "$CROSS_SPECIES_SUMMARY_FILE"; echo "" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "HUMAN TO MOUSE CONSERVATION RATES (based on mapped regions)" >> "$CROSS_SPECIES_SUMMARY_FILE"; echo "-----------------------------" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human liver mapped regions conserved in mouse liver peaks: $HL_C_ML_C - ${HL_C_ML_PCT_S2}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human liver mapped regions conserved in mouse pancreas peaks: $HL_C_MP_C - ${HL_C_MP_PCT_S2}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human pancreas mapped regions conserved in mouse pancreas peaks: $HP_C_MP_C - ${HP_C_MP_PCT_S2}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Human pancreas mapped regions conserved in mouse liver peaks: $HP_C_ML_C - ${HP_C_ML_PCT_S2}%" >> "$CROSS_SPECIES_SUMMARY_FILE"; echo "" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "MOUSE TO HUMAN CONSERVATION STATISTICS" >> "$CROSS_SPECIES_SUMMARY_FILE"; echo "-----------------------------" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Total mouse liver regions: $ML_COUNT_S2" >> "$CROSS_SPECIES_SUMMARY_FILE"; echo "Total mouse pancreas regions: $MP_COUNT_S2" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Mouse liver regions mapped to human: $ML_MAP_C_S2 - ${ML_LIFT_PCT_S2}%" >> "$CROSS_SPECIES_SUMMARY_FILE"; echo "Mouse pancreas regions mapped to human: $MP_MAP_C_S2 - ${MP_LIFT_PCT_S2}%" >> "$CROSS_SPECIES_SUMMARY_FILE"; echo "" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "MOUSE TO HUMAN CONSERVATION RATES (based on mapped regions)" >> "$CROSS_SPECIES_SUMMARY_FILE"; echo "-----------------------------" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Mouse liver mapped regions conserved in human liver peaks: $ML_C_HL_C - ${ML_C_HL_PCT_S2}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Mouse liver mapped regions conserved in human pancreas peaks: $ML_C_HP_C - ${ML_C_HP_PCT_S2}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Mouse pancreas mapped regions conserved in human pancreas peaks: $MP_C_HP_C - ${MP_C_HP_PCT_S2}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "Mouse pancreas mapped regions conserved in human liver peaks: $MP_C_HL_C - ${MP_C_HL_PCT_S2}%" >> "$CROSS_SPECIES_SUMMARY_FILE"
echo "" >> "$CROSS_SPECIES_SUMMARY_FILE"

echo "  Cross-species comparison complete. Results in $CROSS_SPECIES_DIR. Summary in $CROSS_SPECIES_SUMMARY_FILE"
echo "===== Completed Step 2 =====" 