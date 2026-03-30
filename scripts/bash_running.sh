#!/bin/bash

# Get input file paths from command line arguments
INPUT_FILE1=$1
INPUT_FILE2=$2
INPUT_FILE3=$3
INPUT_FILE4=$4
INPUT_FILE5=$5
INPUT_FILE6=$6
INPUT_FILE7=$7
INPUT_FILE8=$8

# Check if all required input files were provided
if [ -z "$INPUT_FILE1" ] || [ -z "$INPUT_FILE2" ] || [ -z "$INPUT_FILE3" ] || [ -z "$INPUT_FILE4" ] || [ -z "$INPUT_FILE5" ] || [ -z "$INPUT_FILE6" ] || [ -z "$INPUT_FILE7" ] || [ -z "$INPUT_FILE8" ]; then
    echo "Error: Not all input files were provided!"
    echo "Usage: ./bash_running.sh <human_liver_peaks> <human_pancreas_peaks> <mouse_liver_peaks> <mouse_pancreas_peaks> <human_liver_to_mouse_halper> <human_pancreas_to_mouse_halper> <mouse_liver_to_human_halper> <mouse_pancreas_to_human_halper>"
    exit 1
fi

# Check if input files exist
for file in "$INPUT_FILE1" "$INPUT_FILE2" "$INPUT_FILE3" "$INPUT_FILE4" "$INPUT_FILE5" "$INPUT_FILE6" "$INPUT_FILE7" "$INPUT_FILE8"; do
    if [ ! -f "$file" ]; then
        echo "Error: Input file $file does not exist"
        exit 1
    fi
done

echo "Running cross-species/tissue mapping"
bash cross_species.sh "$INPUT_FILE1" "$INPUT_FILE2" "$INPUT_FILE3" "$INPUT_FILE4" "$INPUT_FILE5" "$INPUT_FILE6" "$INPUT_FILE7" "$INPUT_FILE8"

# Check if the upstream script executed successfully
if [ $? -ne 0 ]; then
    echo "Error in cross_species script! Pipeline stopped."
    exit 1
fi

# Now define the output files from the upstream script that will be inputs to the R scripts
HUMAN_SHARED="./results/tissue_comparison/human_shared_between_liver_and_pancreas.bed"
HUMAN_LIVER_SPECIFIC="./results/tissue_comparison/human_liver_specific.bed"
HUMAN_PANCREAS_SPECIFIC="./results/tissue_comparison/human_pancreas_specific.bed"
MOUSE_SHARED="./results/tissue_comparison/mouse_shared_between_liver_and_pancreas.bed"
MOUSE_LIVER_SPECIFIC="./results/tissue_comparison/mouse_liver_specific.bed"
MOUSE_PANCREAS_SPECIFIC="./results/tissue_comparison/mouse_pancreas_specific.bed"

HUMAN_LIVER_CONSERVED="./results/cross_species/human_liver_conserved_in_mouse_liver.bed"
HUMAN_LIVER_NOT_CONSERVED="./results/cross_species/human_liver_not_conserved_in_mouse.bed"
MOUSE_LIVER_NOT_CONSERVED="./results/cross_species/mouse_liver_not_conserved_in_human.bed"
HUMAN_PANCREAS_CONSERVED="./results/cross_species/human_pancreas_conserved_in_mouse_pancreas.bed"
HUMAN_PANCREAS_NOT_CONSERVED="./results/cross_species/human_pancreas_not_conserved_in_mouse.bed"
MOUSE_PANCREAS_NOT_CONSERVED="./results/cross_species/mouse_pancreas_not_conserved_in_human.bed"

# Create output directories for R script results
mkdir -p cross-tissue_results cross-species_results

# Step 2: Run cross-tissue R script
echo "Running cross-tissue script"
Rscript CHIPseeker_Cross_tissue.R "$HUMAN_SHARED" "$HUMAN_LIVER_SPECIFIC" "$HUMAN_PANCREAS_SPECIFIC" "$MOUSE_SHARED" "$MOUSE_LIVER_SPECIFIC" "$MOUSE_PANCREAS_SPECIFIC"

# Check if the R script executed successfully
if [ $? -ne 0 ]; then
    echo "Error in cross-tissue script! Pipeline stopped."
    exit 1
fi

# Step 3: Run cross-species R script
echo "Running cross-species script"
Rscript CHIPseeker_Cross_species.R "$HUMAN_LIVER_CONSERVED" "$HUMAN_LIVER_NOT_CONSERVED" "$MOUSE_LIVER_NOT_CONSERVED" "$HUMAN_PANCREAS_CONSERVED" "$HUMAN_PANCREAS_NOT_CONSERVED" "$MOUSE_PANCREAS_NOT_CONSERVED"

# Check if the R script executed successfully
if [ $? -ne 0 ]; then
    echo "Error in cross-species script! Pipeline stopped."
    exit 1
fi

echo "Pipeline completed successfully!"
echo "Results can be found in:"
echo "  - cross-tissue_results (tissue comparison results)"
echo "  - cross-species_results (cross-species conservation results)"