#!/usr/bin/env bash
#SBATCH -p RM-shared
#SBATCH -t 08:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=8
#SBATCH -J halper_liftover
#SBATCH -o /ocean/projects/bio230007p/peerzade/logs/halper_liftover_%j.out
#SBATCH -e /ocean/projects/bio230007p/peerzade/logs/halper_liftover_%j.err

set -euo pipefail

# ─── CONFIGURATION ──────────────────────────────────────────────────────────────
# Project root (for outputs)
PROJECT_ROOT="/ocean/projects/bio230007p/peerzade"

# Shared data root (where HumanAtac/ MouseAtac/ Alignments/ live)
DATA_ROOT="/ocean/projects/bio230007p/ikaplow"

# Cactus HAL alignment
ALIGNMENT_FILE="$DATA_ROOT/Alignments/10plusway-master.hal"

# HALPER repo (holds halper_map_peak_orthologs.sh)
HALPER_SCRIPT_DIR="/jet/home/peerzade/repos/halLiftover-postprocessing"

# Conda environment with halTools + Python deps
CONDA_ENV="halper_env"

# Output directory for mapped peaks
OUTPUT_DIR="$PROJECT_ROOT/mapped_peaks"

# Tissues to process
tissues=( "Liver" "Pancreas" )
# ────────────────────────────────────────────────────────────────────────────────

# Create required directories
mkdir -p "$OUTPUT_DIR" "$PROJECT_ROOT/logs"

# 1) Activate Conda
source /jet/home/peerzade/anaconda3/etc/profile.d/conda.sh
conda activate "$CONDA_ENV"

# 2) Ensure halTools & HALPER scripts are on PATH / PYTHONPATH
export PATH="$HOME/repos/hal/bin:$PATH"
export PYTHONPATH="$HALPER_SCRIPT_DIR:$PYTHONPATH"

# Loop over Human to Mouse and Mouse to Human
for SRC in Human Mouse; do
  TGT=$([[ $SRC == "Human" ]] && echo Mouse || echo Human)

  for TISSUE in "${tissues[@]}"; do
    TAG="${SRC,,}_${TISSUE,,}"
    echo "[$(date +%T)] Mapping ${SRC}→${TGT} for ${TISSUE}"

    # a) Decompress peaks to /tmp
    RAW_GZ="$DATA_ROOT/${SRC}Atac/${TISSUE}/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"
    TEMP_BED="/tmp/${TAG}_peaks.bed"
    zcat "$RAW_GZ" > "$TEMP_BED"

    # b) Run the HALPER wrapper
    bash "$HALPER_SCRIPT_DIR/halper_map_peak_orthologs.sh" \
      -b "$TEMP_BED" \
      -o "$OUTPUT_DIR" \
      -s "$SRC" \
      -t "$TGT" \
      -c "$ALIGNMENT_FILE" \
      -n "$TAG"

    # c) Extract chr/start/end into a simple BED
    MAPPED_NP_GZ="$OUTPUT_DIR/${TAG}.${SRC}To${TGT}.HALPER.narrowPeak.gz"
    zcat "$MAPPED_NP_GZ" | cut -f1-3 \
      > "$OUTPUT_DIR/${TAG}_${TGT,,}_mapped.bed"

    echo "[$(date +%T)] → Wrote:"
    echo " $MAPPED_NP_GZ"
    echo " $OUTPUT_DIR/${TAG}_${TGT,,}_mapped.bed"
  done
done

# Final summary listing
echo
echo "Mapping complete! Contents of $OUTPUT_DIR:"
echo
echo "Raw HALPER output (.HALPER.narrowPeak.gz):"
ls -1 "$OUTPUT_DIR"/*.HALPER.narrowPeak.gz || echo "  (none found!)"
echo
echo "3-col BEDs (_mapped.bed):"
ls -1 "$OUTPUT_DIR"/*_mapped.bed || echo "  (none found!)"
