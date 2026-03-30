#!/bin/bash
#SBATCH --job-name=pipeline
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=64GB

# This script orchestrates the modular pipeline by calling individual step scripts.

set -e

# Load necessary modules
module load bedtools/2.30.0 # Ensure specific version if needed
module load R/4.1.2 # Assuming R is needed for CHIPseeker scripts
module load MEME-suite/5.4.1
module load samtools/1.15 # Specify version if necessary
module load anaconda3/2022.10 # Specify version if necessary
source activate hal # make sure halper environment is named `hal`

# ===== CONFIGURATION: MODIFY THIS ACCORDING TO YOUR PATHS =====
BASE_DIR="/path/to/your/base/directory"
# REFER TO README
# HALPER CONFIG
# Ensure HALPER scripts are executable and in the PATH or specify full path below
HALPER_DIR="/path/to/halLiftover-postprocessing"
export PATH=/path/to/hal/bin:${PATH}
export PYTHONPATH=/path/to/halLiftover-postprocessing:${PYTHONPATH}


# GTF files
HUMAN_GTF="$BASE_DIR/HumanGenomeInfo/gencode.v47.annotation.gff3.gz"
MOUSE_GTF="$BASE_DIR/MouseGenomeInfo/gencode.vM10.annotation.gff3.gz"

# Raw peak files (narrowPeak.gz)
HUMAN_LIVER_PEAKS_GZ="$BASE_DIR/HumanAtac/Liver/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"
HUMAN_PANCREAS_PEAKS_GZ="$BASE_DIR/HumanAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"
MOUSE_LIVER_PEAKS_GZ="$BASE_DIR/MouseAtac/Liver/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"
MOUSE_PANCREAS_PEAKS_GZ="$BASE_DIR/MouseAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"

# Genome FASTA files
HUMAN_GENOME_FA="$BASE_DIR/HumanGenomeInfo/hg38.fa"
MOUSE_GENOME_FA="$BASE_DIR/MouseGenomeInfo/mm10.fa"

# Motif databases
MOTIF_DB_HUMAN="$BASE_DIR/CIS-BP_2.00/Homo_sapiens.meme"
MOTIF_DB_MOUSE="$BASE_DIR/CIS-BP_2.00/Mus_musculus.meme"

# Multi-species alignment
CACTUS_ALIGNMENT="$BASE_DIR/Alignments/10plusway-master.hal"

# Path to the directory containing the modular scripts (assuming they are in the same directory as main.sh)
SCRIPT_DIR=$(dirname "$0")

# ===== CONFIGURATION END =====

# Define top-level output directory with timestamp
OUTPUT_DIR="$BASE_DIR/results_pipeline_$(date +%F_%H-%M)"

# Define expected sub-directory paths based on OUTPUT_DIR (individual scripts will create them)
LOG_DIR="$OUTPUT_DIR/logs"
MAPPED_DIR="$OUTPUT_DIR/mapped_peaks"
STEP2_OUTPUT_DIR="$OUTPUT_DIR/initial_comparisons"
CHIPSEEKER_OUTPUT_DIR="$OUTPUT_DIR/chipseeker_results"
ANALYSIS_OUTPUT_DIR="$OUTPUT_DIR/classification_and_comparisons"
SEQUENCE_DIR="$OUTPUT_DIR/sequences"
MEME_RESULTS="$OUTPUT_DIR/meme_chip_results"

# Create the main output directory and log directory
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# Redirect stdout and stderr for the rest of the script
exec > >(tee -a "$LOG_DIR/main_pipeline_run.log") 2> >(tee -a "$LOG_DIR/main_pipeline_run.err" >&2)

echo "Main pipeline started at: $(date)"
echo "Output directory: $OUTPUT_DIR"
echo "Script directory: $SCRIPT_DIR"

# Index genomes if index files don't exist (can stay in main script or moved to a setup script)
if [ ! -f "${HUMAN_GENOME_FA}.fai" ]; then
  echo "Indexing Human Genome..."
  samtools faidx "$HUMAN_GENOME_FA"
fi
if [ ! -f "${MOUSE_GENOME_FA}.fai" ]; then
  echo "Indexing Mouse Genome..."
  samtools faidx "$MOUSE_GENOME_FA"
fi

# ===== STEP 1: Map with HALPER =====
echo "===== Running Step 1: HALPER Liftover ====="
bash "$SCRIPT_DIR/halper.sh" \
    "$OUTPUT_DIR" \
    "$HALPER_DIR" \
    "$CACTUS_ALIGNMENT" \
    "$HUMAN_LIVER_PEAKS_GZ" \
    "$HUMAN_PANCREAS_PEAKS_GZ" \
    "$MOUSE_LIVER_PEAKS_GZ" \
    "$MOUSE_PANCREAS_PEAKS_GZ"
if [ $? -ne 0 ]; then echo "Error in Step 1" >&2; exit 1; fi
echo "===== Completed Step 1 ====="

# Define expected output paths from Step 1 (needed for subsequent steps)
MAPPED_HL_TO_MM_PEAKS_BED="$MAPPED_DIR/human_liver_Mouse_mapped.bed"
MAPPED_HP_TO_MM_PEAKS_BED="$MAPPED_DIR/human_pancreas_Mouse_mapped.bed"
MAPPED_ML_TO_HG_PEAKS_BED="$MAPPED_DIR/mouse_liver_Human_mapped.bed"
MAPPED_MP_TO_HG_PEAKS_BED="$MAPPED_DIR/mouse_pancreas_Human_mapped.bed"
HUMAN_LIVER_TO_MOUSE_HALPER_NP="$MAPPED_DIR/human_liver.HumanToMouse.HALPER.narrowPeak.gz"
HUMAN_PANCREAS_TO_MOUSE_HALPER_NP="$MAPPED_DIR/human_pancreas.HumanToMouse.HALPER.narrowPeak.gz"
MOUSE_LIVER_TO_HUMAN_HALPER_NP="$MAPPED_DIR/mouse_liver.MouseToHuman.HALPER.narrowPeak.gz"
MOUSE_PANCREAS_TO_HUMAN_HALPER_NP="$MAPPED_DIR/mouse_pancreas.MouseToHuman.HALPER.narrowPeak.gz"

# ===== STEP 2: Initial Cross-species and Tissue Comparison =====
echo "===== Running Step 2: Initial Bedtools Comparisons ====="
bash "$SCRIPT_DIR/compare.sh" \
    "$OUTPUT_DIR" \
    "$HUMAN_LIVER_PEAKS_GZ" \
    "$HUMAN_PANCREAS_PEAKS_GZ" \
    "$MOUSE_LIVER_PEAKS_GZ" \
    "$MOUSE_PANCREAS_PEAKS_GZ" \
    "$HUMAN_LIVER_TO_MOUSE_HALPER_NP" \
    "$HUMAN_PANCREAS_TO_MOUSE_HALPER_NP" \
    "$MOUSE_LIVER_TO_HUMAN_HALPER_NP" \
    "$MOUSE_PANCREAS_TO_HUMAN_HALPER_NP"
if [ $? -ne 0 ]; then echo "Error in Step 2" >&2; exit 1; fi
echo "===== Completed Step 2 ====="

# Define expected output paths from Step 2 (needed for Step 3)
# Tissue Comparison Outputs
HUMAN_SHARED_STEP2_BED="$STEP2_OUTPUT_DIR/tissue_comparison/human_shared_between_liver_and_pancreas.bed"
HUMAN_LIVER_SPECIFIC_STEP2_BED="$STEP2_OUTPUT_DIR/tissue_comparison/human_liver_specific.bed"
HUMAN_PANCREAS_SPECIFIC_STEP2_BED="$STEP2_OUTPUT_DIR/tissue_comparison/human_pancreas_specific.bed"
MOUSE_SHARED_STEP2_BED="$STEP2_OUTPUT_DIR/tissue_comparison/mouse_shared_between_liver_and_pancreas.bed"
MOUSE_LIVER_SPECIFIC_STEP2_BED="$STEP2_OUTPUT_DIR/tissue_comparison/mouse_liver_specific.bed"
MOUSE_PANCREAS_SPECIFIC_STEP2_BED="$STEP2_OUTPUT_DIR/tissue_comparison/mouse_pancreas_specific.bed"
# Cross-Species Comparison Outputs (Mapped coordinates in target species)
HL_CONS_ML_BED="$STEP2_OUTPUT_DIR/cross_species/human_liver_conserved_in_mouse_liver.bed"
HL_NOT_CONS_S2_BED="$STEP2_OUTPUT_DIR/cross_species/human_liver_not_conserved_in_mouse.bed"
ML_NOT_CONS_S2_BED="$STEP2_OUTPUT_DIR/cross_species/mouse_liver_not_conserved_in_human.bed"
HP_CONS_MP_BED="$STEP2_OUTPUT_DIR/cross_species/human_pancreas_conserved_in_mouse_pancreas.bed"
HP_NOT_CONS_S2_BED="$STEP2_OUTPUT_DIR/cross_species/human_pancreas_not_conserved_in_mouse.bed"
MP_NOT_CONS_S2_BED="$STEP2_OUTPUT_DIR/cross_species/mouse_pancreas_not_conserved_in_human.bed"

# ===== STEP 3: Run CHIPseeker R Scripts =====
echo "===== Running Step 3: CHIPseeker Analysis ====="
# Ensure R scripts are findable (e.g., in SCRIPT_DIR or PATH)
# Assuming R scripts are in the same directory as main.sh
bash "$SCRIPT_DIR/chipseeker.sh" \
    "$OUTPUT_DIR" \
    "$HUMAN_SHARED_STEP2_BED" \
    "$HUMAN_LIVER_SPECIFIC_STEP2_BED" \
    "$HUMAN_PANCREAS_SPECIFIC_STEP2_BED" \
    "$MOUSE_SHARED_STEP2_BED" \
    "$MOUSE_LIVER_SPECIFIC_STEP2_BED" \
    "$MOUSE_PANCREAS_SPECIFIC_STEP2_BED" \
    "$HL_CONS_ML_BED" \
    "$HL_NOT_CONS_S2_BED" \
    "$ML_NOT_CONS_S2_BED" \
    "$HP_CONS_MP_BED" \
    "$HP_NOT_CONS_S2_BED" \
    "$MP_NOT_CONS_S2_BED"
if [ $? -ne 0 ]; then echo "Error in Step 3" >&2; exit 1; fi
echo "===== Completed Step 3 ====="

# ===== STEP 4: Classify Peaks, Compare Subsets & Calculate Stats =====
echo "===== Running Step 4: Classification and Comparisons ====="
bash "$SCRIPT_DIR/pro_enh.sh" \
    "$OUTPUT_DIR" \
    "$HUMAN_GTF" \
    "$MOUSE_GTF" \
    "$HUMAN_LIVER_PEAKS_GZ" \
    "$HUMAN_PANCREAS_PEAKS_GZ" \
    "$MOUSE_LIVER_PEAKS_GZ" \
    "$MOUSE_PANCREAS_PEAKS_GZ" \
    "$MAPPED_HL_TO_MM_PEAKS_BED" \
    "$MAPPED_HP_TO_MM_PEAKS_BED" \
    "$MAPPED_ML_TO_HG_PEAKS_BED" \
    "$MAPPED_MP_TO_HG_PEAKS_BED"
if [ $? -ne 0 ]; then echo "Error in Step 4" >&2; exit 1; fi
echo "===== Completed Step 4 ====="

# ===== STEP 5: Run MEME-ChIP Analysis =====
echo "===== Running Step 5: MEME-ChIP ====="
bash "$SCRIPT_DIR/memechip.sh" \
    "$OUTPUT_DIR" \
    "$HUMAN_LIVER_PEAKS_GZ" \
    "$HUMAN_PANCREAS_PEAKS_GZ" \
    "$MOUSE_LIVER_PEAKS_GZ" \
    "$MOUSE_PANCREAS_PEAKS_GZ" \
    "$HUMAN_GENOME_FA" \
    "$MOUSE_GENOME_FA" \
    "$MOTIF_DB_HUMAN" \
    "$MOTIF_DB_MOUSE" \
    "$MAPPED_HL_TO_MM_PEAKS_BED" \
    "$MAPPED_HP_TO_MM_PEAKS_BED" \
    "$MAPPED_ML_TO_HG_PEAKS_BED" \
    "$MAPPED_MP_TO_HG_PEAKS_BED"
if [ $? -ne 0 ]; then echo "Error in Step 5" >&2; exit 1; fi
echo "===== Completed Step 5 ====="


echo "===== Main pipeline completed successfully at: $(date) ====="