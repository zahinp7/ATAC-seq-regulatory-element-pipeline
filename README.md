# Cross-Species Regulatory Element Analysis Pipeline

**Team:** Andres Chousal, Zahin Peerzade, Siddharth Sabata, Yinuo Yang

**Course:** Bioinformatics Data Practicum

**Tissues:** Liver, Pancreas (Human \& Mouse)

---

## Introduction

This pipeline performs a comparative analysis of transcriptional regulatory elements (promoters and enhancers) between liver and pancreas tissues in human and mouse. It maps open chromatin regions across species, identifies conserved and tissue-specific elements, classifies peaks as promoters or enhancers, and performs motif discovery and functional enrichment.

We aim to answer three primary questions with these results:

**Question 1**: Is transcriptional regulatory element activity more conserved across tissues or species?

**Question 2**: To what extent does the transcriptional regulatory code differ between enhancers and promoters?

**Question 3**: To what extent are the biological processes upregulated in tissue conserved across species?

(*Questions directly quoted from project description document*)

---
## Workflow

![](https://github.com/achousal/Bioinformatics_03713/blob/main/pipeline_final.png)
1. **Data Quality Control**
ATAC-seq (Assay for Transposase-Accessible Chromatin using Sequencing) quality reports (provided) from human and mouse liver, pancreas, and ovary manually analyzed to determine which datasets to use for further analysis. The ovarian dataset consisted of short read lengths for both human and mouse, and was thrown out. This step is essential for ensuring informative results. 
2. **Cross-Species Mapping**
Open chromatin regions identified from ATAC-seq are mapped between human and mouse genomes using [HALPER](https://github.com/pfenninglab/halLiftover-postprocessing) (HAL Liftover Post-processing for Epigenomic Regions) with Cactus whole-genome alignments. This step is crucial for identifying orthologous regulatory elements between species, allowing direct comparison of regulatory activity at corresponding genomic locations. 
3. **Conservation and Specificity Analysis (Question 1)**
This step identifies regulatory elements based on their conservation patterns: elements conserved between species for the same tissue, elements shared across tissues within a species, and elements specific to a tissue or species. By quantifying these different categories, we can directly address whether regulatory element activity is more conserved across tissues or species. 
4. **Functional Annotation (Question 3)** 
Regulatory elements identified in previous steps are annotated and analyzed for Gene Ontology (GO) term enrichment with [ChIPseeker](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html). This analysis connects regulatory elements to their potential target genes and biological functions, revealing which biological processes are regulated by conserved versus species-specific elements. The functional analysis helps determine whether similar biological processes are regulated by conserved elements across species, despite potential differences in the specific regulatory elements.
5. **Enhancer/Promoter Classification (Question 2)** 
Open chromatin regions are classified as promoters (within 2kb upstream and 200bp downstream of Tanscription Start Site, or TSS) or enhancers (all other regions) using [BEDTools](https://bedtools.readthedocs.io/) and genome annotations. This classification is essential for understanding how conservation patterns differ between these two types of regulatory elements. Promoters, which are proximal to genes, may be under different evolutionary constraints than enhancers, which can act at a distance and may evolve more rapidly.
6. **Motif Discovery (Question 2)** 
Sequences from classified regulatory elements are extracted and analyzed using [MEME-ChIP](https://meme-suite.org/meme/) for de novo motif discovery. This step identifies enriched transcription factor binding motifs in different categories of regulatory elements, revealing how the transcriptional regulatory code differs between tissues and species. By comparing motifs between enhancers and promoters across species and tissues, we can understand how the language of transcription factors has evolved.

---
## Repository Setup

- `pipeline/` contains an automated, end to end script of this pipeline. It has NOT been tested due to technical difficulties with the Pittsburgh Supercomputing Center. 
    - Scripts in this directory explained later in this document.
- `scripts/` contains prototypes, and original scripts used to generate data. While doing this project, we used smaller scripts to run each part instead of one end to end pipeline. These individual scripts were adapted into the end to end pipeline. 
    - `bash_running.sh`: Bash script for running ChIPseeker.
    - `cross_species.sh`: Slurm script for cross-species/tissue analysis with BEDTools.
    - `enhancer.sh`: Slurm script for enhancer/promoter classification and statistics.
    - `halper_liftover.sh`: Slurm script for HALPER.
    - `map_halper.sh`: Slurm script for HALPER.
    - `run_full_pipeline.sh`: Slurm script prototype of full pipeline. Final pipeline took this implementation, but split into six files. 
    - `CHIPseeker/`: ChIPseeker R scripts, along with associated reports. 
    - `meme_chip/`: MEME-ChIP results. 

---
## Setup and Configuration (`pipeline/main.sh`)

The main pipeline execution is controlled by `main.sh`. This script handles the overall setup, configuration, and sequential execution of the individual step scripts. Before running the pipeline, ensure the following setup and configuration steps are correctly addressed:

**Execution Environment & Dependencies:**

**Requirements**
- [HALPER](https://github.com/pfenninglab/halLiftover-postprocessing)
    - Follow [these](https://github.com/pfenninglab/halLiftover-postprocessing/blob/master/hal_install_instructions.md) directions exactly
- [BEDTools](https://bedtools.readthedocs.io/)
- [MEME Suite](https://meme-suite.org/)
- [ChIPseeker (R package)](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
- [Anaconda3](https://www.anaconda.com/)
- [R](https://www.r-project.org/)
- [samtools](https://www.htslib.org/)
- HPC with SLURM

Install dependencies using conda. Closely follow installation instructions for installing each required tool. 

**Required Input Files** 
- **ATAC-seq peak files:**
    - 4 total: human liver, human pancreas, mouse liver, mouse pancreas
    - `idr.optimal_peak.narrowPeak.gz`
    - Pipeline automatically decrompresses 
- **Cactus (multi-species) alignment:**
    - `alignment.hal`
- **Genome annotations:**
    - 2 total: human, mouse
    - `annotation.gff3.gz`
    - Pipeline automatically decrompresses 
- **Reference genomes:**
    - 2 total: human, mouse
    - `ref_genome.fa`
- **Motif database:**
    - 2 total: human, mouse
    - `motif_db.meme`

All paths to input files and key directories must be correctly set within the `CONFIGURATION` section of `main.sh`:

*   `BASE_DIR`: The top-level base directory. Input files are typically expected relative to this directory, and the main output directory will be created here.
*   `HALPER_DIR`: Path to the directory containing HALPER helper scripts, specifically `halper_map_peak_orthologs.sh` used in Step 1.
    * For the `HALPER CONFIG` step, make sure to alter the following accordingly in the configuration section of `main.sh`:

    ```
    export PATH=[repos dir]/hal/bin:${PATH}
    export PYTHONPATH=[repos dir]/halLiftover-postprocessing:${PYTHONPATH}
    ```

*   `HUMAN_GTF`, `MOUSE_GTF`: Full paths to the gzipped gene annotation files (GFF3 format) for human (hg38) and mouse (mm10).
*   `*_PEAKS_GZ` (e.g., `HUMAN_LIVER_PEAKS_GZ`): Full paths to the original gzipped narrowPeak files containing ATAC-seq peaks for each species and tissue.
*   `*_GENOME_FA` (e.g., `HUMAN_GENOME_FA`): Full paths to the reference genome FASTA files for human and mouse.
*   `MOTIF_DB_HUMAN`, `MOTIF_DB_MOUSE`: Full paths to the motif databases in MEME format for human and mouse, used by MEME-ChIP in Step 5.
*   `CACTUS_ALIGNMENT`: Full path to the multi-species whole-genome alignment file in HAL format, used by HALPER in Step 1.

**Script Locations:**

*   `SCRIPT_DIR`: This variable is automatically set to the directory where `main.sh` is located. The script assumes that all the step scripts (`halper.sh`, `compare.sh`, `chipseeker.sh`, `pro_enh.sh`, `memechip.sh`) and the R scripts (`CHIPseeker_Cross_tissue.R`, `CHIPseeker_Cross_species.R`) reside in this same directory. If they are located elsewhere, you will need to adjust the paths in the `bash` and `Rscript` calls within `main.sh` and `chipseeker.sh`.

**Output Directory:**

*   `OUTPUT_DIR`: An output directory named `results_pipeline_YYYY-MM-DD_HH-MM` (timestamped) will be created under `BASE_DIR`. All results, logs, and intermediate files from the pipeline steps will be organized within this directory.
*   `LOG_DIR`: A `logs` subdirectory within `OUTPUT_DIR` will store the main stdout (`main_pipeline_run.log`) and stderr (`main_pipeline_run.err`) for the entire pipeline run. Individual step scripts might create their own logs if configured to do so (currently they don't).

**5. Genome Indexing:**

*   The `main.sh` script checks for the existence of FASTA index files (`.fai`) for the human and mouse genomes. If an index file is missing, it automatically generates it using `samtools faidx`. Ensure you have write permissions in the directory containing the genome FASTA files if indexing is required.

**6. Running the Pipeline:**

*   Submit the `main.sh` script to your scheduler (e.g., `sbatch main.sh`) or run it directly in an interactive session after ensuring all dependencies and configurations are correct.
*   The script executes each step sequentially. If any step fails, the pipeline will exit with an error message.

The main pipeline is orchestrated by `main.sh`, which calls the following scripts in sequence.

---

## Step 1: `halper.sh` - HALPER Liftover

*   **Purpose**: Maps peaks between species (Human <-> Mouse) using the HALPER tool. Generates mapped coordinates in BED3 format while retaining the original detailed HALPER output in narrowPeak format.
*   **Usage**:
    ```bash
    ./halper.sh <output_dir> <halper_dir> <cactus_alignment> <human_liver_peaks_gz> <human_pancreas_peaks_gz> <mouse_liver_peaks_gz> <mouse_pancreas_peaks_gz>
    ```
*   **Inputs**:
    *   `Arg 1: <output_dir>`: Base directory for all pipeline results.
    *   `Arg 2: <halper_dir>`: Directory containing HALPER helper scripts (e.g., `halper_map_peak_orthologs.sh`).
    *   `Arg 3: <cactus_alignment>`: Path to the HAL alignment file (e.g., `*.hal`).
    *   `Arg 4: <human_liver_peaks_gz>`: Gzipped narrowPeak - Original Human Liver peaks.
    *   `Arg 5: <human_pancreas_peaks_gz>`: Gzipped narrowPeak - Original Human Pancreas peaks.
    *   `Arg 6: <mouse_liver_peaks_gz>`: Gzipped narrowPeak - Original Mouse Liver peaks.
    *   `Arg 7: <mouse_pancreas_peaks_gz>`: Gzipped narrowPeak - Original Mouse Pancreas peaks.
*   **Outputs** (within `<output_dir>/mapped_peaks/`):
    *   `human_liver_Mouse_mapped.bed`: BED3 - Human Liver peaks mapped to Mouse coordinates.
    *   `human_pancreas_Mouse_mapped.bed`: BED3 - Human Pancreas peaks mapped to Mouse coordinates.
    *   `mouse_liver_Human_mapped.bed`: BED3 - Mouse Liver peaks mapped to Human coordinates.
    *   `mouse_pancreas_Human_mapped.bed`: BED3 - Mouse Pancreas peaks mapped to Human coordinates.
    *   `human_liver.HumanToMouse.HALPER.narrowPeak.gz`: Gzipped narrowPeak - Raw HALPER output (Human Liver -> Mouse).
    *   `human_pancreas.HumanToMouse.HALPER.narrowPeak.gz`: Gzipped narrowPeak - Raw HALPER output (Human Pancreas -> Mouse).
    *   `mouse_liver.MouseToHuman.HALPER.narrowPeak.gz`: Gzipped narrowPeak - Raw HALPER output (Mouse Liver -> Human).
    *   `mouse_pancreas.MouseToHuman.HALPER.narrowPeak.gz`: Gzipped narrowPeak - Raw HALPER output (Mouse Pancreas -> Human).

---

## Step 2: `compare.sh` - Initial Bedtools Comparisons

*   **Purpose**: Performs initial tissue-level (shared/specific) and cross-species conservation comparisons using `bedtools`. It processes the raw HALPER narrowPeak output from Step 1 and compares regions against the original peak sets.
*   **Usage**:
    ```bash
    ./compare.sh <output_dir> <human_liver_peaks_gz> <human_pancreas_peaks_gz> <mouse_liver_peaks_gz> <mouse_pancreas_peaks_gz> <human_liver_to_mouse_halper_np> <human_pancreas_to_mouse_halper_np> <mouse_liver_to_human_halper_np> <mouse_pancreas_to_human_halper_np>
    ```
*   **Inputs**:
    *   `Arg 1: <output_dir>`: Base directory for pipeline results.
    *   `Arg 2: <human_liver_peaks_gz>`: Gzipped narrowPeak - Original Human Liver peaks.
    *   `Arg 3: <human_pancreas_peaks_gz>`: Gzipped narrowPeak - Original Human Pancreas peaks.
    *   `Arg 4: <mouse_liver_peaks_gz>`: Gzipped narrowPeak - Original Mouse Liver peaks.
    *   `Arg 5: <mouse_pancreas_peaks_gz>`: Gzipped narrowPeak - Original Mouse Pancreas peaks.
    *   `Arg 6: <human_liver_to_mouse_halper_np>`: Gzipped narrowPeak - HALPER output (Human Liver -> Mouse) from Step 1.
    *   `Arg 7: <human_pancreas_to_mouse_halper_np>`: Gzipped narrowPeak - HALPER output (Human Pancreas -> Mouse) from Step 1.
    *   `Arg 8: <mouse_liver_to_human_halper_np>`: Gzipped narrowPeak - HALPER output (Mouse Liver -> Human) from Step 1.
    *   `Arg 9: <mouse_pancreas_to_human_halper_np>`: Gzipped narrowPeak - HALPER output (Mouse Pancreas -> Human) from Step 1.
*   **Outputs** (within `<output_dir>/initial_comparisons/`):
    *   **Tissue Comparison** (`tissue_comparison/`):
        *   BED files: `human_shared_between_liver_and_pancreas.bed`, `human_liver_specific.bed`, `human_pancreas_specific.bed`, `mouse_shared_between_liver_and_pancreas.bed`, `mouse_liver_specific.bed`, `mouse_pancreas_specific.bed`.
        *   Summary: `step2_tissue_comparison_summary.txt`.
    *   **Cross-Species Comparison** (`cross_species/`):
        *   Processed HALPER BEDs (`processed_files/`): `human_liver_to_mouse.bed`, etc.
        *   Conserved BEDs: `human_liver_conserved_in_mouse_liver.bed`, `human_liver_conserved_in_mouse_pancreas.bed`, etc. (8 files total).
        *   Non-Conserved BEDs: `human_liver_not_conserved_in_mouse.bed`, `human_pancreas_not_conserved_in_mouse.bed`, etc. (4 files total).
        *   Summary: `step2_cross_species_conservation_summary.txt`.

---

## Step 3: `chipseeker.sh` - CHIPseeker Annotation

*   **Purpose**: Runs CHIPseeker R scripts to perform genomic feature annotation on the peak subsets generated in Step 2.
*   **Usage**:
    ```bash
    ./chipseeker.sh <output_dir> <human_shared> <human_liver_specific> <human_pancreas_specific> <mouse_shared> <mouse_liver_specific> <mouse_pancreas_specific> <human_liver_conserved> <human_liver_not_conserved> <mouse_liver_not_conserved> <human_pancreas_conserved> <human_pancreas_not_conserved> <mouse_pancreas_not_conserved>
    ```
*   **Inputs**:
    *   `Arg 1: <output_dir>`: Base directory for pipeline results.
    *   `Arg 2-7`: BED files for tissue shared/specific peaks (Human/Mouse) from Step 2.
    *   `Arg 8-13`: BED files for cross-species conserved/not-conserved peaks (Human/Mouse) from Step 2.
    *   Dependencies: `CHIPseeker_Cross_tissue.R`, `CHIPseeker_Cross_species.R` (R scripts assumed to be accessible).
*   **Outputs** (within `<output_dir>/chipseeker_results/`):
    *   `cross_tissue_results/`: Directory containing annotation plots and tables from `CHIPseeker_Cross_tissue.R`.
    *   `cross_species_results/`: Directory containing annotation plots and tables from `CHIPseeker_Cross_species.R`.

---

## Step 4: `pro_enh.sh` - Promoter/Enhancer Classification & Comparison

*   **Purpose**: Defines promoter regions, classifies original peaks into promoter or enhancer subsets, performs detailed tissue and species comparisons on these classified subsets, and calculates summary statistics.
*   **Usage**:
    ```bash
    ./pro_enh.sh <output_dir> <human_gtf> <mouse_gtf> <human_liver_peaks_gz> <human_pancreas_peaks_gz> <mouse_liver_peaks_gz> <mouse_pancreas_peaks_gz> <mapped_hl_to_mm_bed> <mapped_hp_to_mm_bed> <mapped_ml_to_hg_bed> <mapped_mp_to_hg_bed>
    ```
*   **Inputs**:
    *   `Arg 1: <output_dir>`: Base directory for pipeline results.
    *   `Arg 2: <human_gtf>`: Gzipped GFF3 - Human gene annotation.
    *   `Arg 3: <mouse_gtf>`: Gzipped GFF3 - Mouse gene annotation.
    *   `Arg 4-7`: Gzipped narrowPeak files - Original peak sets (Human/Mouse, Liver/Pancreas).
    *   `Arg 8-11`: BED files - Mapped peak coordinates (from Step 1, e.g., `human_liver_Mouse_mapped.bed`).
*   **Outputs** (within `<output_dir>/classification_and_comparisons/`):
    *   **Classification** (`classified_peaks/`):
        *   Promoter Definitions: `human_promoters_definition.bed`, `mouse_promoters_definition.bed`.
        *   Classified Peaks (Human): `human/human_liver_promoters.bed`, `human/human_liver_enhancers.bed`, etc. (4 files).
        *   Classified Peaks (Mouse): `mouse/mouse_liver_promoters.bed`, `mouse/mouse_liver_enhancers.bed`, etc. (4 files).
    *   **Tissue Comparison** (`region_comparisons/tissue_comparison/`):
        *   Shared/Specific Enhancers: `human_enhancers_tissue_shared.bed`, `human_liver_enhancers_tissue_specific.bed`, etc. (6 files).
        *   Shared/Specific Promoters: `human_promoters_tissue_shared.bed`, `human_liver_promoters_tissue_specific.bed`, etc. (6 files).
    *   **Species Comparison** (`region_comparisons/species_comparison/`):
        *   Shared/Specific Enhancers vs Mapped: `human_liver_enhancers_shared_with_mapped_mouse_peaks_human_coords.bed`, `human_liver_enhancers_specific_vs_mapped_mouse_peaks_human_coords.bed`, etc. (8 files).
        *   Shared/Specific Promoters vs Mapped: `human_liver_promoters_shared_with_mapped_mouse_peaks_human_coords.bed`, `human_liver_promoters_specific_vs_mapped_mouse_peaks_human_coords.bed`, etc. (8 files).
    *   **Summary Statistics**: `classification_and_comparison_summary.txt`.

---

## Step 5: `memechip.sh` - MEME-ChIP Motif Analysis

*   **Purpose**: Performs motif discovery using MEME-ChIP on various peak subsets generated in Step 4 (including full peak sets, classified promoters/enhancers, and tissue/species specific/shared enhancers).
*   **Usage**:
    ```bash
    ./memechip.sh <output_dir> <human_liver_peaks_gz> <human_pancreas_peaks_gz> <mouse_liver_peaks_gz> <mouse_pancreas_peaks_gz> <human_genome_fa> <mouse_genome_fa> <motif_db_human> <motif_db_mouse> <mapped_hl_to_mm_bed> <mapped_hp_to_mm_bed> <mapped_ml_to_hg_bed> <mapped_mp_to_hg_bed>
    ```
*   **Inputs**:
    *   `Arg 1: <output_dir>`: Base directory for pipeline results.
    *   `Arg 2-5`: Gzipped narrowPeak - Original peak sets (used for 'full' subset analysis).
    *   `Arg 6: <human_genome_fa>`: FASTA - Human reference genome.
    *   `Arg 7: <mouse_genome_fa>`: FASTA - Mouse reference genome.
    *   `Arg 8: <motif_db_human>`: MEME format - Human motif database.
    *   `Arg 9: <motif_db_mouse>`: MEME format - Mouse motif database.
    *   `Arg 10-13`: BED - Mapped peaks from Step 1 (used only to check if species comparison was skipped in Step 4).
    *   *Implicit Dependencies*: Requires classified/comparison BED files generated by `pro_enh.sh` (Step 4) to be present in `<output_dir>/classification_and_comparisons/`.
*   **Outputs**:
    *   `<output_dir>/classification_and_comparisons/full_peak_beds/*.bed`: BED3 versions of original peaks (created by this script).
    *   `<output_dir>/sequences/<subset_name>.fa`: FASTA sequences for each analyzed subset.
    *   `<output_dir>/meme_chip_results/<subset_name>/`: MEME-ChIP output directory for each subset, containing standard results like `meme-chip.html`, `meme_out/`, etc.
    *   `<output_dir>/meme_chip_results/<subset_name>/meme-chip.log`: Log file for each MEME-ChIP run.

--- 
## Cite Us

If you use this pipeline in your research, please cite:

Andres Chousal, Zahin Peerzade, Siddharth Sabata, Yinuo Yang. *RegulatoryElementAnalysisPipeline: A pipeline for Cross-Species Regulatory Element Analysis of Human and Mouse ATAC data*. GitHub repository. Available at: https://github.com/BioinformaticsDataPracticum2025/Bioinformatics_03713 (Accessed 29 Apr 2025).


---

## References

- HALPER: Xiaoyu Zhang, Irene Kaplow, Morgan Wirthlin, Tyler Park, Andreas Pfenning. HALPER facilitates the identification of regulatory element orthologs across species. *Bioinformatics*, Volume 36, Issue 15, 1 August 2020, Pages 4339-4340. 
- MEME Suite: Bailey TL, Johnson J, Grant CE, Noble WS. The MEME Suite. *Nucleic Acids Res*. 2015;43(W1):W39-W49. doi:10.1093/nar/gkv416 
- ChIPseeker: Yu G, Wang LG, He QY. ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. *Bioinformatics*. 2015;31(14):2382-2383. doi:10.1093/bioinformatics/btv145 
- BEDTools: Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*. 2010;26(6):841-842. doi:10.1093/bioinformatics/btq033 

---

## AI Usage
- GenAI assistance used to create README
- GenAI used in all code. Original code was written by humans, but AI was used to optimize, comment, and bugfix. The final scripts in this repository have been modified with AI. 
