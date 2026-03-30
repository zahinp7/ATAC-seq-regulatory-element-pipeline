#!/bin/bash
# halper.sh - Map peaks between species using HALPER.
#             Generates mapped coordinates (BED3) and keeps original HALPER output (narrowPeak.gz).
#
# Usage: ./halper.sh <output_dir> <halper_dir> <cactus_alignment> <human_liver_peaks_gz> <human_pancreas_peaks_gz> <mouse_liver_peaks_gz> <mouse_pancreas_peaks_gz>
#
# Inputs:
#   Arg 1: <output_dir>              - Base directory for pipeline results.
#   Arg 2: <halper_dir>              - Directory containing HALPER helper scripts (e.g., halper_map_peak_orthologs.sh).
#   Arg 3: <cactus_alignment>        - Path to the HAL alignment file (e.g., *.hal).
#   Arg 4: <human_liver_peaks_gz>    - Gzipped narrowPeak: Original Human Liver peaks.
#   Arg 5: <human_pancreas_peaks_gz> - Gzipped narrowPeak: Original Human Pancreas peaks.
#   Arg 6: <mouse_liver_peaks_gz>    - Gzipped narrowPeak: Original Mouse Liver peaks.
#   Arg 7: <mouse_pancreas_peaks_gz> - Gzipped narrowPeak: Original Mouse Pancreas peaks.
#
# Outputs (within <output_dir>/mapped_peaks/ directory):
#   File: human_liver_Mouse_mapped.bed                - BED3: Human Liver peaks mapped to Mouse coordinates.
#   File: human_pancreas_Mouse_mapped.bed             - BED3: Human Pancreas peaks mapped to Mouse coordinates.
#   File: mouse_liver_Human_mapped.bed                - BED3: Mouse Liver peaks mapped to Human coordinates.
#   File: mouse_pancreas_Human_mapped.bed             - BED3: Mouse Pancreas peaks mapped to Human coordinates.
#   File: human_liver.HumanToMouse.HALPER.narrowPeak.gz   - Gzipped narrowPeak: Raw HALPER output for Human Liver -> Mouse.
#   File: human_pancreas.HumanToMouse.HALPER.narrowPeak.gz - Gzipped narrowPeak: Raw HALPER output for Human Pancreas -> Mouse.
#   File: mouse_liver.MouseToHuman.HALPER.narrowPeak.gz   - Gzipped narrowPeak: Raw HALPER output for Mouse Liver -> Human.
#   File: mouse_pancreas.MouseToHuman.HALPER.narrowPeak.gz - Gzipped narrowPeak: Raw HALPER output for Mouse Pancreas -> Human.

set -e

# Parse command-line arguments
if [ $# -lt 7 ]; then
    echo "Usage: $0 <output_dir> <halper_dir> <cactus_alignment> <human_liver_peaks_gz> <human_pancreas_peaks_gz> <mouse_liver_peaks_gz> <mouse_pancreas_peaks_gz>"
    exit 1
fi

OUTPUT_DIR=$1
HALPER_DIR=$2
CACTUS_ALIGNMENT=$3
HUMAN_LIVER_PEAKS=$4
HUMAN_PANCREAS_PEAKS=$5
MOUSE_LIVER_PEAKS=$6
MOUSE_PANCREAS_PEAKS=$7

# Create output directory
MAPPED_DIR="$OUTPUT_DIR/mapped_peaks"
LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$LOG_DIR" "$MAPPED_DIR"

# Function to run HALPER and create a simple BED3 output
function run_halper() {
  local species_source=$1; local species_target=$2; local input_peaks_gz=$3; local name=$4
  local output_mapped_bed="$MAPPED_DIR/${name}_${species_target}_mapped.bed" # BED3 Mapped peaks
  local halper_output_np_path="$MAPPED_DIR/${name}.${species_source}To${species_target}.HALPER.narrowPeak.gz" # Keep track of original narrowPeak output
  local tmp_bed="/tmp/${name}_peaks_$$.bed"

  # Check if BED3 output already exists
  if [ -s "$output_mapped_bed" ]; then
      echo "  Skipping HALPER for $name ($species_source to $species_target): Output $output_mapped_bed already exists."
      # Ensure the narrowPeak output also exists if skipping
      if [ ! -f "$halper_output_np_path" ]; then
          echo "Error: BED3 output $output_mapped_bed exists, but original HALPER narrowPeak $halper_output_np_path is missing. Please check HALPER step." >&2
          exit 1
      fi
      return
  fi

  echo "  Running HALPER for $name ($species_source to $species_target)...";
  zcat "$input_peaks_gz" > "$tmp_bed"; if [ ! -s "$tmp_bed" ]; then echo "Error: Failed to create or empty temporary BED file $tmp_bed from $input_peaks_gz" >&2; rm -f "$tmp_bed"; exit 1; fi
  # Run HALPER script (assuming it produces the .narrowPeak.gz file)
  bash "$HALPER_DIR/halper_map_peak_orthologs.sh" -b "$tmp_bed" -o "$MAPPED_DIR" -s "$species_source" -t "$species_target" -c "$CACTUS_ALIGNMENT" -n "$name";
  # Check if narrowPeak output was created
  if [ ! -f "$halper_output_np_path" ]; then echo "Error: HALPER output file $halper_output_np_path not found for $name ($species_source to $species_target)." >&2; rm -f "$tmp_bed"; exit 1; fi
  # Process narrowPeak to BED3
  zcat "$halper_output_np_path" | cut -f1-3 > "$output_mapped_bed"; if [ ! -s "$output_mapped_bed" ]; then echo "Warning: Created empty mapped BED file $output_mapped_bed for $name ($species_source to $species_target)." >&2; fi
  rm "$tmp_bed"; echo "  Finished HALPER for $name ($species_source to $species_target)."
}

echo "===== STEP 1: Running HALPER liftover ====="

run_halper Human Mouse "$HUMAN_LIVER_PEAKS" human_liver
run_halper Human Mouse "$HUMAN_PANCREAS_PEAKS" human_pancreas
run_halper Mouse Human "$MOUSE_LIVER_PEAKS" mouse_liver
run_halper Mouse Human "$MOUSE_PANCREAS_PEAKS" mouse_pancreas

echo "Step 1 complete. Mapped peaks are in $MAPPED_DIR." 