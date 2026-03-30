#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 04:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=9
#SBATCH -J halper_mapping
#SBATCH -o halper_mapping_%j.out
#SBATCH -e halper_mapping_%j.err

# Create logs directory if it doesn't exist
mkdir -p logs

# Set paths
HALPER_DIR="/path/to/halLiftover-postprocessing"
OUTPUT_DIR="$PROJECT/mapped_peaks"
CACTUS_ALIGNMENT="/path/to/your/base/directory/Alignments/10plusway-master.hal"

# Create output directory
mkdir -p $OUTPUT_DIR

# Activate conda environment with HALPER dependencies
module load anaconda3
conda init
source activate hal

export PATH=/path/to/hal/bin:${PATH}
export PYTHONPATH=/path/to/halLiftover-postprocessing:${PYTHONPATH}

# Human Liver
zcat /path/to/your/base/directory/HumanAtac/Liver/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz > /tmp/human_liver_peaks.bed
bash $HALPER_DIR/halper_map_peak_orthologs.sh \
  -b /tmp/human_liver_peaks.bed \
  -o $OUTPUT_DIR \
  -s Human \
  -t Mouse \
  -c $CACTUS_ALIGNMENT \
  -n human_liver
zcat $OUTPUT_DIR/human_liver.HumanToMouse.HALPER.narrowPeak.gz | cut -f1-3 > $OUTPUT_DIR/human_liver_mouse_mapped.bed

# Human Pancreas
zcat /path/to/your/base/directory/HumanAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz > /tmp/human_pancreas_peaks.bed
bash $HALPER_DIR/halper_map_peak_orthologs.sh \
  -b /tmp/human_pancreas_peaks.bed \
  -o $OUTPUT_DIR \
  -s Human \
  -t Mouse \
  -c $CACTUS_ALIGNMENT \
  -n human_pancreas
zcat $OUTPUT_DIR/human_pancreas.HumanToMouse.HALPER.narrowPeak.gz | cut -f1-3 > $OUTPUT_DIR/human_pancreas_mouse_mapped.bed

# Mouse Liver
zcat /path/to/your/base/directory/MouseAtac/Liver/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz > /tmp/mouse_liver_peaks.bed
bash $HALPER_DIR/halper_map_peak_orthologs.sh \
  -b /tmp/mouse_liver_peaks.bed \
  -o $OUTPUT_DIR \
  -s Mouse \
  -t Human \
  -c $CACTUS_ALIGNMENT \
  -n mouse_liver
zcat $OUTPUT_DIR/mouse_liver.MouseToHuman.HALPER.narrowPeak.gz | cut -f1-3 > $OUTPUT_DIR/mouse_liver_human_mapped.bed

# Mouse Pancreas
zcat /path/to/your/base/directory/MouseAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz > /tmp/mouse_pancreas_peaks.bed
bash $HALPER_DIR/halper_map_peak_orthologs.sh \
  -b /tmp/mouse_pancreas_peaks.bed \
  -o $OUTPUT_DIR \
  -s Mouse \
  -t Human \
  -c $CACTUS_ALIGNMENT \
  -n mouse_pancreas
zcat $OUTPUT_DIR/mouse_pancreas.MouseToHuman.HALPER.narrowPeak.gz | cut -f1-3 > $OUTPUT_DIR/mouse_pancreas_human_mapped.bed

echo "All HALPER mappings and GREAT-ready BED files created."