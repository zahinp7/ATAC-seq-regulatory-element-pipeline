#!/bin/bash
# chipseeker.sh - Run CHIPseeker R scripts for genomic annotation analysis
#
# Usage: ./chipseeker.sh <output_dir> <human_shared> <human_liver_specific> <human_pancreas_specific> <mouse_shared> <mouse_liver_specific> <mouse_pancreas_specific> <human_liver_conserved> <human_liver_not_conserved> <mouse_liver_not_conserved> <human_pancreas_conserved> <human_pancreas_not_conserved> <mouse_pancreas_not_conserved>
#
# Inputs:
#   Arg 1:  <output_dir>                     - Base directory for pipeline results.
#   Arg 2:  <human_shared>                   - BED: Human peaks shared between liver & pancreas (Step 2 output).
#   Arg 3:  <human_liver_specific>           - BED: Human liver-specific peaks (Step 2 output).
#   Arg 4:  <human_pancreas_specific>        - BED: Human pancreas-specific peaks (Step 2 output).
#   Arg 5:  <mouse_shared>                   - BED: Mouse peaks shared between liver & pancreas (Step 2 output).
#   Arg 6:  <mouse_liver_specific>           - BED: Mouse liver-specific peaks (Step 2 output).
#   Arg 7:  <mouse_pancreas_specific>        - BED: Mouse pancreas-specific peaks (Step 2 output).
#   Arg 8:  <human_liver_conserved>          - BED: Mapped human liver peaks conserved in mouse liver (Step 2 output).
#   Arg 9:  <human_liver_not_conserved>      - BED: Mapped human liver peaks not conserved in any mouse tissue (Step 2 output).
#   Arg 10: <mouse_liver_not_conserved>      - BED: Mapped mouse liver peaks not conserved in any human tissue (Step 2 output).
#   Arg 11: <human_pancreas_conserved>       - BED: Mapped human pancreas peaks conserved in mouse pancreas (Step 2 output).
#   Arg 12: <human_pancreas_not_conserved>   - BED: Mapped human pancreas peaks not conserved in any mouse tissue (Step 2 output).
#   Arg 13: <mouse_pancreas_not_conserved>   - BED: Mapped mouse pancreas peaks not conserved in any human tissue (Step 2 output).
#   Dependency: CHIPseeker_Cross_tissue.R    - R script performing cross-tissue annotation.
#   Dependency: CHIPseeker_Cross_species.R   - R script performing cross-species annotation.
#
# Outputs:
#   - Directory: <output_dir>/chipseeker_results/cross_tissue_results/   - Contains annotation plots and tables from CHIPseeker_Cross_tissue.R.
#   - Directory: <output_dir>/chipseeker_results/cross_species_results/ - Contains annotation plots and tables from CHIPseeker_Cross_species.R.

set -e

# Parse command-line arguments
if [ $# -lt 13 ]; then
    echo "Usage: $0 <output_dir> <human_shared> <human_liver_specific> <human_pancreas_specific> <mouse_shared> <mouse_liver_specific> <mouse_pancreas_specific> <human_liver_conserved> <human_liver_not_conserved> <mouse_liver_not_conserved> <human_pancreas_conserved> <human_pancreas_not_conserved> <mouse_pancreas_not_conserved>"
    exit 1
fi

OUTPUT_DIR=$1
HUMAN_SHARED=$2
HUMAN_LIVER_SPECIFIC=$3
HUMAN_PANCREAS_SPECIFIC=$4
MOUSE_SHARED=$5
MOUSE_LIVER_SPECIFIC=$6
MOUSE_PANCREAS_SPECIFIC=$7
HUMAN_LIVER_CONSERVED=$8
HUMAN_LIVER_NOT_CONSERVED=$9
MOUSE_LIVER_NOT_CONSERVED=${10}
HUMAN_PANCREAS_CONSERVED=${11}
HUMAN_PANCREAS_NOT_CONSERVED=${12}
MOUSE_PANCREAS_NOT_CONSERVED=${13}

# Create output directory
CHIPSEEKER_OUTPUT_DIR="$OUTPUT_DIR/chipseeker_results"
mkdir -p "$CHIPSEEKER_OUTPUT_DIR" "$CHIPSEEKER_OUTPUT_DIR/cross_tissue_results" "$CHIPSEEKER_OUTPUT_DIR/cross_species_results"

echo "===== STEP 3: Running CHIPseeker R scripts ====="

# Verify BED files for R scripts exist
echo "  Verifying input BED files for R scripts..."
all_r_inputs_found=true
for bed_file in "$HUMAN_SHARED" "$HUMAN_LIVER_SPECIFIC" "$HUMAN_PANCREAS_SPECIFIC" \
                "$MOUSE_SHARED" "$MOUSE_LIVER_SPECIFIC" "$MOUSE_PANCREAS_SPECIFIC" \
                "$HUMAN_LIVER_CONSERVED" "$HUMAN_LIVER_NOT_CONSERVED" "$MOUSE_LIVER_NOT_CONSERVED" \
                "$HUMAN_PANCREAS_CONSERVED" "$HUMAN_PANCREAS_NOT_CONSERVED" "$MOUSE_PANCREAS_NOT_CONSERVED"; do
    if [ ! -s "$bed_file" ]; then
        echo "Error: Required BED file $bed_file for R scripts not found or empty." >&2
        all_r_inputs_found=false
    fi
done
if ! $all_r_inputs_found ; then
    exit 1
fi
echo "  All input BED files for R scripts found."

# Cross-tissue R analysis
echo "  Running CHIPseeker_Cross_tissue.R..."
Rscript CHIPseeker_Cross_tissue.R "$HUMAN_SHARED" "$HUMAN_LIVER_SPECIFIC" "$HUMAN_PANCREAS_SPECIFIC" \
                                  "$MOUSE_SHARED" "$MOUSE_LIVER_SPECIFIC" "$MOUSE_PANCREAS_SPECIFIC" \
                                  --output_dir "$CHIPSEEKER_OUTPUT_DIR/cross_tissue_results"

if [ $? -ne 0 ]; then
  echo "Error in CHIPseeker_Cross_tissue.R" >&2
  exit 1
fi

# Cross-species R analysis
echo "  Running CHIPseeker_Cross_species.R..."
Rscript CHIPseeker_Cross_species.R "$HUMAN_LIVER_CONSERVED" "$HUMAN_LIVER_NOT_CONSERVED" "$MOUSE_LIVER_NOT_CONSERVED" \
                                   "$HUMAN_PANCREAS_CONSERVED" "$HUMAN_PANCREAS_NOT_CONSERVED" "$MOUSE_PANCREAS_NOT_CONSERVED" \
                                   --output_dir "$CHIPSEEKER_OUTPUT_DIR/cross_species_results"

if [ $? -ne 0 ]; then
  echo "Error in CHIPseeker_Cross_species.R" >&2
  exit 1
fi

echo "  CHIPseeker analysis complete. Results in: $CHIPSEEKER_OUTPUT_DIR"
echo "===== Completed Step 3 =====" 