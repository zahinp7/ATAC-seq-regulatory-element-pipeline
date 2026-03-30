#!/bin/bash
#SBATCH --job-name=full_pipeline
#SBATCH --output=logs/full_pipeline_%j.out
#SBATCH --error=logs/full_pipeline_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=64GB

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
export PATH=/path/to/hal/bin:${PATH}
export PYTHONPATH=/path/to/halLiftover-postprocessing:${PYTHONPATH}
HALPER_DIR="/path/to/halLiftover-postprocessing"

# GTF files
HUMAN_GTF="$BASE_DIR/HumanGenomeInfo/gencode.v47.annotation.gff3.gz"
MOUSE_GTF="$BASE_DIR/MouseGenomeInfo/gencode.vM10.annotation.gff3.gz"

# Raw peak files
HUMAN_LIVER_PEAKS="$BASE_DIR/HumanAtac/Liver/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"
HUMAN_PANCREAS_PEAKS="$BASE_DIR/HumanAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"
MOUSE_LIVER_PEAKS="$BASE_DIR/MouseAtac/Liver/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"
MOUSE_PANCREAS_PEAKS="$BASE_DIR/MouseAtac/Pancreas/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz"

# Genome FASTA files
HUMAN_GENOME="$BASE_DIR/HumanGenomeInfo/hg38.fa"
MOUSE_GENOME="$BASE_DIR/MouseGenomeInfo/mm10.fa"

# Motif databases
MOTIF_DB_HUMAN="$BASE_DIR/CIS-BP_2.00/Homo_sapiens.meme"
MOTIF_DB_MOUSE="$BASE_DIR/CIS-BP_2.00/Mus_musculus.meme"

# Multi-species alignment
CACTUS_ALIGNMENT="$BASE_DIR/Alignments/10plusway-master.hal"

# ===== CONFIGURATION END =====

#top-level output directories
OUTPUT_DIR="$BASE_DIR/results_pipeline_$(date +%F_%H-%M)" # Add timestamp to avoid overwriting
MAPPED_DIR="$OUTPUT_DIR/mapped_peaks" # Directory for HALPER mapped BED files
STEP2_OUTPUT_DIR="$OUTPUT_DIR/initial_comparisons" # Directory for output of the new Step 2
CHIPSEEKER_OUTPUT_DIR="$OUTPUT_DIR/chipseeker_results" # Dedicated dir for CHIPseeker outputs
ANALYSIS_OUTPUT_DIR="$OUTPUT_DIR/classification_and_comparisons" # Dedicated dir for classification/comparison step output
SEQUENCE_DIR="$OUTPUT_DIR/sequences" # Directory for FASTA sequences for MEME
MEME_RESULTS="$OUTPUT_DIR/meme_chip_results" # Directory for MEME-ChIP outputs
LOG_DIR="$OUTPUT_DIR/logs"
# Create all necessary directories upfront
mkdir -p "$LOG_DIR" "$OUTPUT_DIR" "$MAPPED_DIR" "$STEP2_OUTPUT_DIR" "$CHIPSEEKER_OUTPUT_DIR" \
         "$ANALYSIS_OUTPUT_DIR" "$SEQUENCE_DIR" "$MEME_RESULTS"
mkdir -p "$ANALYSIS_OUTPUT_DIR/classified_peaks/human" "$ANALYSIS_OUTPUT_DIR/classified_peaks/mouse"
mkdir -p "$ANALYSIS_OUTPUT_DIR/region_comparisons/tissue_comparison" "$ANALYSIS_OUTPUT_DIR/region_comparisons/species_comparison"


# Redirect stdout and stderr for the rest of the script
exec > >(tee -a "$LOG_DIR/pipeline_run.log") 2> >(tee -a "$LOG_DIR/pipeline_run.err" >&2)

echo "Pipeline started at: $(date)"
echo "Output directory: $OUTPUT_DIR"

# Index genomes if index files don't exist
if [ ! -f "${HUMAN_GENOME}.fai" ]; then
  echo "Indexing Human Genome..."
samtools faidx "$HUMAN_GENOME"
fi
if [ ! -f "${MOUSE_GENOME}.fai" ]; then
  echo "Indexing Mouse Genome..."
samtools faidx "$MOUSE_GENOME"
fi

# ===== STEP 1: Map with HALPER =====
echo "===== STEP 1: Running HALPER liftover ====="

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

run_halper Human Mouse "$HUMAN_LIVER_PEAKS" human_liver
run_halper Human Mouse "$HUMAN_PANCREAS_PEAKS" human_pancreas
run_halper Mouse Human "$MOUSE_LIVER_PEAKS" mouse_liver
run_halper Mouse Human "$MOUSE_PANCREAS_PEAKS" mouse_pancreas

# Define paths to the mapped BED files created by HALPER run (Used in Step 4)
MAPPED_HL_TO_MM_PEAKS="$MAPPED_DIR/human_liver_Mouse_mapped.bed"
MAPPED_HP_TO_MM_PEAKS="$MAPPED_DIR/human_pancreas_Mouse_mapped.bed"
MAPPED_ML_TO_HG_PEAKS="$MAPPED_DIR/mouse_liver_Human_mapped.bed"
MAPPED_MP_TO_HG_PEAKS="$MAPPED_DIR/mouse_pancreas_Human_mapped.bed"

# Define paths to the original HALPER narrowPeak outputs (Used in Step 2)
HUMAN_LIVER_TO_MOUSE_HALPER="$MAPPED_DIR/human_liver.HumanToMouse.HALPER.narrowPeak.gz"
HUMAN_PANCREAS_TO_MOUSE_HALPER="$MAPPED_DIR/human_pancreas.HumanToMouse.HALPER.narrowPeak.gz"
MOUSE_LIVER_TO_HUMAN_HALPER="$MAPPED_DIR/mouse_liver.MouseToHuman.HALPER.narrowPeak.gz"
MOUSE_PANCREAS_TO_HUMAN_HALPER="$MAPPED_DIR/mouse_pancreas.MouseToHuman.HALPER.narrowPeak.gz"

# ===== STEP 2: Initial Cross-species and Tissue Comparison (bedtools based) =====

echo "===== STEP 2: Running Initial Bedtools Comparisons ====="

# Define output subdirectories for this step
CROSS_SPECIES_DIR="$STEP2_OUTPUT_DIR/cross_species"
TISSUE_COMP_DIR="$STEP2_OUTPUT_DIR/tissue_comparison"
PROCESSED_DIR="$CROSS_SPECIES_DIR/processed_files" # For processed HALPER outputs
mkdir -p "$CROSS_SPECIES_DIR" "$TISSUE_COMP_DIR" "$PROCESSED_DIR"

# Log input files used in this step
echo "  Inputs for Step 2:"
echo "    Human liver peaks: $HUMAN_LIVER_PEAKS"
echo "    Human pancreas peaks: $HUMAN_PANCREAS_PEAKS"
echo "    Mouse liver peaks: $MOUSE_LIVER_PEAKS"
echo "    Mouse pancreas peaks: $MOUSE_PANCREAS_PEAKS"
echo "    Human liver to mouse HALPER NP: $HUMAN_LIVER_TO_MOUSE_HALPER"
echo "    Human pancreas to mouse HALPER NP: $HUMAN_PANCREAS_TO_MOUSE_HALPER"
echo "    Mouse liver to human HALPER NP: $MOUSE_LIVER_TO_HUMAN_HALPER"
echo "    Mouse pancreas to human HALPER NP: $MOUSE_PANCREAS_TO_HUMAN_HALPER"

# Check if *input* files for this step exist (Original Peaks and HALPER narrowPeaks)
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

# --- PART 2.1: WITHIN-SPECIES TISSUE COMPARISON (using original peaks) ---
echo "  Starting Tissue Comparison Analysis (using original peaks)..."

# Convert narrowPeak files to BED format (more specific BED6-like format)
echo "  Converting narrowPeak files to BED format..."
HUMAN_LIVER_BED="$TISSUE_COMP_DIR/human_liver.bed"
HUMAN_PANCREAS_BED="$TISSUE_COMP_DIR/human_pancreas.bed"
MOUSE_LIVER_BED="$TISSUE_COMP_DIR/mouse_liver.bed"
MOUSE_PANCREAS_BED="$TISSUE_COMP_DIR/mouse_pancreas.bed"
zcat "$HUMAN_LIVER_PEAKS" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,".",$7}' > "$HUMAN_LIVER_BED" 
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
TISSUE_SUMMARY_FILE="$STEP2_OUTPUT_DIR/step2_tissue_comparison_summary.txt" # Specific summary file
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

# ===== STEP 3: Run CHIPseeker R Scripts =====
echo "===== STEP 3: Running CHIPseeker R scripts ====="

# Define BED inputs for R analysis, pointing to the output of *this* Step 2
HUMAN_SHARED="$HUMAN_SHARED_STEP2"
HUMAN_LIVER_SPECIFIC="$HUMAN_LIVER_SPECIFIC_STEP2"
HUMAN_PANCREAS_SPECIFIC="$HUMAN_PANCREAS_SPECIFIC_STEP2"
MOUSE_SHARED="$MOUSE_SHARED_STEP2"
MOUSE_LIVER_SPECIFIC="$MOUSE_LIVER_SPECIFIC_STEP2"
MOUSE_PANCREAS_SPECIFIC="$MOUSE_PANCREAS_SPECIFIC_STEP2"

HUMAN_LIVER_CONSERVED="$HL_CONS_ML" # Human liver peaks whose mapped region overlaps mouse liver peaks
HUMAN_LIVER_NOT_CONSERVED="$HL_NOT_CONS_S2" # Human liver peaks whose mapped region overlaps *no* mouse peaks
MOUSE_LIVER_NOT_CONSERVED="$ML_NOT_CONS_S2" # Mouse liver peaks whose mapped region overlaps *no* human peaks
HUMAN_PANCREAS_CONSERVED="$HP_CONS_MP" # Human pancreas peaks whose mapped region overlaps mouse pancreas peaks
HUMAN_PANCREAS_NOT_CONSERVED="$HP_NOT_CONS_S2" # Human pancreas peaks whose mapped region overlaps *no* mouse peaks
MOUSE_PANCREAS_NOT_CONSERVED="$MP_NOT_CONS_S2" # Mouse pancreas peaks whose mapped region overlaps *no* human peaks

# Verify BED files for R scripts exist
echo "  Verifying input BED files for R scripts..."
all_r_inputs_found=true
for bed_file in "$HUMAN_SHARED" "$HUMAN_LIVER_SPECIFIC" "$HUMAN_PANCREAS_SPECIFIC" \
                "$MOUSE_SHARED" "$MOUSE_LIVER_SPECIFIC" "$MOUSE_PANCREAS_SPECIFIC" \
                "$HUMAN_LIVER_CONSERVED" "$HUMAN_LIVER_NOT_CONSERVED" "$MOUSE_LIVER_NOT_CONSERVED" \
                "$HUMAN_PANCREAS_CONSERVED" "$HUMAN_PANCREAS_NOT_CONSERVED" "$MOUSE_PANCREAS_NOT_CONSERVED"; do
    if [ ! -s "$bed_file" ]; then
        echo "Error: Required BED file $bed_file for R scripts not found or empty. Check Step 2 outputs in $STEP2_OUTPUT_DIR." >&2
        all_r_inputs_found=false
    fi
done
if ! $all_r_inputs_found ; then
    exit 1
fi
echo "  All input BED files for R scripts found."

# Cross-tissue R analysis - outputting to CHIPSEEKER_OUTPUT_DIR
echo "  Running CHIPseeker_Cross_tissue.R..."
Rscript CHIPseeker_Cross_tissue.R "$HUMAN_SHARED" "$HUMAN_LIVER_SPECIFIC" "$HUMAN_PANCREAS_SPECIFIC" \
                                  "$MOUSE_SHARED" "$MOUSE_LIVER_SPECIFIC" "$MOUSE_PANCREAS_SPECIFIC" \
                                  --output_dir "$CHIPSEEKER_OUTPUT_DIR/cross_tissue_results" 

if [ $? -ne 0 ]; then
  echo "Error in CHIPseeker_Cross_tissue.R" >&2
  exit 1
fi

# Cross-species R analysis - outputting to CHIPSEEKER_OUTPUT_DIR
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


# ===== STEP 4: Classify Peaks, Compare Subsets & Calculate Stats =====
# This step generates the Promoter/Enhancer classifications and comparisons needed for the final questions and MEME-ChIP
echo "===== STEP 4: Classifying Peaks, Comparing Subsets & Calculating Stats ====="

# Define directories
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

# ===== STEP 5: Run MEME-ChIP Analysis =====
echo "===== STEP 5: Running MEME-ChIP motif analysis ====="

FULL_PEAKS_BED_DIR="$ANALYSIS_OUTPUT_DIR/full_peak_beds" # Dir for BED3 version of full peaks
mkdir -p "$FULL_PEAKS_BED_DIR" "$SEQUENCE_DIR" "$MEME_RESULTS" # Ensure sequence/meme dirs exist

echo "  Preparing full peak set BED files (if needed)..."
# Create BED3 files for the full peak sets (category a) only if they don't exist
if [[ ! -s "$FULL_PEAKS_BED_DIR/human_liver_full.bed" ]]; then zcat "$HUMAN_LIVER_PEAKS" | cut -f 1-3 > "$FULL_PEAKS_BED_DIR/human_liver_full.bed"; fi
if [[ ! -s "$FULL_PEAKS_BED_DIR/human_pancreas_full.bed" ]]; then zcat "$HUMAN_PANCREAS_PEAKS" | cut -f 1-3 > "$FULL_PEAKS_BED_DIR/human_pancreas_full.bed"; fi
if [[ ! -s "$FULL_PEAKS_BED_DIR/mouse_liver_full.bed" ]]; then zcat "$MOUSE_LIVER_PEAKS" | cut -f 1-3 > "$FULL_PEAKS_BED_DIR/mouse_liver_full.bed"; fi
if [[ ! -s "$FULL_PEAKS_BED_DIR/mouse_pancreas_full.bed" ]]; then zcat "$MOUSE_PANCREAS_PEAKS" | cut -f 1-3 > "$FULL_PEAKS_BED_DIR/mouse_pancreas_full.bed"; fi

# Define the list of BED files for MEME-ChIP (categories a-g from Step 4 outputs)
declare -A meme_inputs=(
    # a. Full peak sets
    ["human_liver_full"]="$FULL_PEAKS_BED_DIR/human_liver_full.bed"
    ["human_pancreas_full"]="$FULL_PEAKS_BED_DIR/human_pancreas_full.bed"
    ["mouse_liver_full"]="$FULL_PEAKS_BED_DIR/mouse_liver_full.bed"
    ["mouse_pancreas_full"]="$FULL_PEAKS_BED_DIR/mouse_pancreas_full.bed"
    # b. Enhancers
    ["human_liver_enhancers"]="$HUMAN_LIVER_ENHANCERS_BED"
    ["human_pancreas_enhancers"]="$HUMAN_PANCREAS_ENHANCERS_BED"
    ["mouse_liver_enhancers"]="$MOUSE_LIVER_ENHANCERS_BED"
    ["mouse_pancreas_enhancers"]="$MOUSE_PANCREAS_ENHANCERS_BED"
    # c. Promoters
    ["human_liver_promoters"]="$HUMAN_LIVER_PROMOTERS_BED"
    ["human_pancreas_promoters"]="$HUMAN_PANCREAS_PROMOTERS_BED"
    ["mouse_liver_promoters"]="$MOUSE_LIVER_PROMOTERS_BED"
    ["mouse_pancreas_promoters"]="$MOUSE_PANCREAS_PROMOTERS_BED"
    # d. Enhancers shared across tissues
    ["human_enhancers_tissue_shared"]="$HUMAN_ENHANCERS_TISSUE_SHARED"
    ["mouse_enhancers_tissue_shared"]="$MOUSE_ENHANCERS_TISSUE_SHARED"
    # e. Enhancers specific to each tissue
    ["human_liver_enhancers_tissue_specific"]="$HUMAN_LIVER_ENHANCERS_TISSUE_SPECIFIC"
    ["human_pancreas_enhancers_tissue_specific"]="$HUMAN_PANCREAS_ENHANCERS_TISSUE_SPECIFIC"
    ["mouse_liver_enhancers_tissue_specific"]="$MOUSE_LIVER_ENHANCERS_TISSUE_SPECIFIC"
    ["mouse_pancreas_enhancers_tissue_specific"]="$MOUSE_PANCREAS_ENHANCERS_TISSUE_SPECIFIC"
    # f. Enhancers shared across species (vs mapped peaks)
    ["human_liver_enhancers_species_shared_hg"]="$HUMAN_LIVER_ENH_SHARED_MAPPED_MOUSE_HG_COORDS"
    ["human_pancreas_enhancers_species_shared_hg"]="$HUMAN_PANCREAS_ENH_SHARED_MAPPED_MOUSE_HG_COORDS"
    ["mouse_liver_enhancers_species_shared_mm"]="$MOUSE_LIVER_ENH_SHARED_MAPPED_HUMAN_MM_COORDS"
    ["mouse_pancreas_enhancers_species_shared_mm"]="$MOUSE_PANCREAS_ENH_SHARED_MAPPED_HUMAN_MM_COORDS"
    # g. Enhancers specific to each species (vs mapped peaks)
    ["human_liver_enhancers_species_specific_hg"]="$HUMAN_LIVER_ENH_SPECIFIC_MAPPED_MOUSE_HG_COORDS"
    ["human_pancreas_enhancers_species_specific_hg"]="$HUMAN_PANCREAS_ENH_SPECIFIC_MAPPED_MOUSE_HG_COORDS"
    ["mouse_liver_enhancers_species_specific_mm"]="$MOUSE_LIVER_ENH_SPECIFIC_MAPPED_HUMAN_MM_COORDS"
    ["mouse_pancreas_enhancers_species_specific_mm"]="$MOUSE_PANCREAS_ENH_SPECIFIC_MAPPED_HUMAN_MM_COORDS"
)

# Loop through the defined inputs and run MEME-ChIP
for name in "${!meme_inputs[@]}"; do
  BED_FILE="${meme_inputs[$name]}"
  FASTA_FILE="$SEQUENCE_DIR/${name}.fa"
  MEME_OUT_DIR="$MEME_RESULTS/${name}"

  # Determine species and genome/motif DB
  if [[ "$name" == human* ]]; then GENOME="$HUMAN_GENOME"; SPECIES="human"; MOTIF_DB="$MOTIF_DB_HUMAN";
  elif [[ "$name" == mouse* ]]; then GENOME="$MOUSE_GENOME"; SPECIES="mouse"; MOTIF_DB="$MOTIF_DB_MOUSE";
  else echo "    Warning: Cannot determine species for $name based on name. Skipping MEME-ChIP." >&2; continue; fi

  echo "    Processing $name for MEME-ChIP..."

  # Check if BED file exists and is not empty
  # For species comparison files, also check if the step was skipped
  if [[ "$SKIP_SPECIES_STATS" -eq 1 && "$name" == *species_* ]]; then echo "    Skipping $name: Species comparison step was skipped."; continue; fi
  if [[ ! -s "$BED_FILE" ]]; then echo "    Skipping $name: Input BED file $BED_FILE not found or empty."; continue; fi

  # Create FASTA file
  if [[ -s "$FASTA_FILE" ]]; then echo "      Skipping getfasta for $name: FASTA file exists."
  else echo "      Running getfasta for $name..."; bedtools getfasta -fi "$GENOME" -bed "$BED_FILE" -fo "$FASTA_FILE"; if [[ ! -s "$FASTA_FILE" ]]; then echo "    Skipping $name: Failed to create FASTA file or it is empty."; continue; fi; fi

  # Run MEME-ChIP
  if [[ -d "$MEME_OUT_DIR/meme_out" || -f "$MEME_OUT_DIR/meme-chip.html" ]]; then echo "      Skipping MEME-ChIP for $name: Output directory exists."
  else mkdir -p "$MEME_OUT_DIR"; echo "      Running MEME-ChIP for $name..."; meme-chip -oc "$MEME_OUT_DIR" -db "$MOTIF_DB" -meme-nmotifs 5 -meme-maxw 20 -seed 42 "$FASTA_FILE" &> "$MEME_OUT_DIR/meme-chip.log"; if [ $? -ne 0 ]; then echo "Error: MEME-ChIP failed for $name. Check log: $MEME_OUT_DIR/meme-chip.log" >&2; else echo "      Finished MEME-ChIP for $name."; fi; fi
done

echo "===== Completed Step 5 ====="


echo "===== Pipeline completed successfully at: $(date) ====="