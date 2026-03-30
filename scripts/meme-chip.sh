#!/bin/bash
#SBATCH --job-name=full_meme_chip
#SBATCH --output=logs/full_meme_chip_%j.out
#SBATCH --error=logs/full_meme_chip_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=18
#SBATCH --mem=32GB
#SBATCH --partition=RM-shared

# Load modules
module load bedtools
module load MEME-suite/5.4.1
module load samtools  # for faidx

# Directories
BASE_DIR="/ocean/projects/bio230007p/achousal"
CLASSIFIED_DIR="${BASE_DIR}/promoters_enhancers"
GENOME_DIR="${BASE_DIR}/ikaplow"
OUTPUT_BASE="${BASE_DIR}/meme_chip_outputs"
mkdir -p logs "$OUTPUT_BASE"

# Genome FASTA files
HUMAN_GENOME="${GENOME_DIR}/HumanGenomeInfo/hg38.fa"
MOUSE_GENOME="${GENOME_DIR}/MouseGenomeInfo/mm10.fa"

# Ensure FASTA index exists
samtools faidx "$HUMAN_GENOME"
samtools faidx "$MOUSE_GENOME"

# Sample configuration: [species] [tissue] [region_type]
samples=(
  "human liver promoters"
  "human liver enhancers"
  "human pancreas promoters"
  "human pancreas enhancers"
  "mouse liver promoters"
  "mouse liver enhancers"
  "mouse pancreas promoters"
  "mouse pancreas enhancers"
)

# Loop through each sample
for sample in "${samples[@]}"; do
  read -r SPECIES TISSUE REGION <<< "$sample"

  echo "Processing $SPECIES $TISSUE $REGION..."

  # Set paths
  BED_FILE="${CLASSIFIED_DIR}/${SPECIES}/${SPECIES}_${TISSUE}_${REGION}.bed"
  FASTA_FILE="${BASE_DIR}/mapped_peaks/${SPECIES}_${TISSUE}_${REGION}.fa"
  OUTPUT_DIR="${OUTPUT_BASE}/${SPECIES}_${TISSUE}_${REGION}"
  mkdir -p "$OUTPUT_DIR"

  # Select genome and motif DB
  if [[ "$SPECIES" == "human" ]]; then
    GENOME="$HUMAN_GENOME"
    MOTIF_DB="${GENOME_DIR}/CIS-BP_2.00/Homo_sapiens.meme"
  else
    GENOME="$MOUSE_GENOME"
    MOTIF_DB="${GENOME_DIR}/CIS-BP_2.00/Mus_musculus.meme"
  fi

  # Convert BED to FASTA
  bedtools getfasta -fi "$GENOME" -bed "$BED_FILE" -fo "$FASTA_FILE"

  # Check if FASTA file has content
  if [[ ! -s "$FASTA_FILE" ]]; then
    echo "Skipping $SPECIES $TISSUE $REGION: FASTA file is empty or missing."
    continue
  fi

  # Run MEME-ChIP
  meme-chip \
    -oc "$OUTPUT_DIR" \
    -db "$MOTIF_DB" \
    -meme-nmotifs 5 \
    -maxw 20 \
    "$FASTA_FILE"

  echo "Finished $SPECIES $TISSUE $REGION."
done

echo "All MEME-ChIP analyses completed."
