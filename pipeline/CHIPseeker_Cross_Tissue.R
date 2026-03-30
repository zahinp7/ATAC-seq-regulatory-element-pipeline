#############################################################
#This script is for functional annotation of human and mouse#
#conserved regions of liver and pancreas (cross-tissue)######
#############################################################
rm(list = ls())
# Install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChIPseeker")
# For mouse 
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene") # For mouse
BiocManager::install("org.Mm.eg.db") # Mouse gene IDs
# For human 
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")
# Install GO.db
BiocManager::install("GO.db")
# Install clusterProfiler
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
# Additional packages for visualization
BiocManager::install(c("enrichplot", "ggplot2"))

# Load libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Create output directory
dir.create("cross-tissue_results", showWarnings = FALSE)

# Import the peak data (BED file)
args <- commandArgs(trailingOnly = TRUE)
input_file1 <- args[1]
input_file2 <- args[2]
input_file3 <- args[3]
input_file4 <- args[4]
input_file5 <- args[5]
input_file6 <- args[6]

# Human .bed file 
human_conserve_peaks <- readPeakFile(input_file1)
human_liver_specific_peaks <- readPeakFile(input_file2)
human_pancreas_specific_peaks <- readPeakFile(input_file3)

# Mouse .bed file
mouse_conserve_peaks <- readPeakFile(input_file4) 
mouse_liver_specific_peaks <- readPeakFile(input_file5)
mouse_pancreas_specific_peaks <- readPeakFile(input_file6)


# Get a summary of peak annotations
txdb_mouse <- TxDb.Mmusculus.UCSC.mm10.knownGene
txdb_human <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Human 
peakAnno_human_conserve <- annotatePeak(human_conserve_peaks, tssRegion=c(-3000, 3000),
                                        TxDb=txdb_human, annoDb="org.Hs.eg.db")
peakAnno_human_liver <- annotatePeak(human_liver_specific_peaks, tssRegion=c(-3000, 3000),
                                     TxDb=txdb_human, annoDb="org.Hs.eg.db")
peakAnno_human_pancreas <- annotatePeak(human_pancreas_specific_peaks, tssRegion=c(-3000, 3000),
                                        TxDb=txdb_human, annoDb="org.Hs.eg.db")
# Mouse 
peakAnno_mouse_conserve <- annotatePeak(mouse_conserve_peaks, tssRegion=c(-3000, 3000),
                                        TxDb=txdb_mouse, annoDb="org.Mm.eg.db")
peakAnno_mouse_liver <- annotatePeak(mouse_liver_specific_peaks, tssRegion=c(-3000, 3000),
                                     TxDb=txdb_mouse, annoDb="org.Mm.eg.db")
peakAnno_mouse_pancreas <- annotatePeak(mouse_pancreas_specific_peaks, tssRegion=c(-3000, 3000),
                                        TxDb=txdb_mouse, annoDb="org.Mm.eg.db")

# View annotation summary
summary(peakAnno_human_conserve)
summary(peakAnno_human_liver)
summary(peakAnno_human_pancreas)
summary(peakAnno_mouse_conserve)
summary(peakAnno_mouse_liver)
summary(peakAnno_mouse_pancreas)

# Convert to data frame for easier manipulation
peakAnnoDF_human_conserve <- as.data.frame(peakAnno_human_conserve)
peakAnnoDF_human_liver <- as.data.frame(peakAnno_human_liver)
peakAnnoDF_human_pancreas <- as.data.frame(peakAnno_human_pancreas)
peakAnnoDF_mouse_conserve <- as.data.frame(peakAnno_mouse_conserve)
peakAnnoDF_mouse_liver <- as.data.frame(peakAnno_mouse_liver)
peakAnnoDF_mouse_pancreas <- as.data.frame(peakAnno_mouse_pancreas)

# Create a list of all annotated peaks for combined visualization
all_peaks <- list(
  "Human_Conserved" = peakAnno_human_conserve,
  "Human_Liver" = peakAnno_human_liver,
  "Human_Pancreas" = peakAnno_human_pancreas,
  "Mouse_Conserved" = peakAnno_mouse_conserve,
  "Mouse_Liver" = peakAnno_mouse_liver,
  "Mouse_Pancreas" = peakAnno_mouse_pancreas
)


# Combined bar plot of genomic feature distribution
pdf("cross-tissue_results/combined_anno_bar.pdf", width=10, height=6)
plotAnnoBar(all_peaks)
dev.off()

# Combined distribution of peaks relative to TSS
pdf("cross-tissue_results/combined_distToTSS.pdf", width=10, height=6)
plotDistToTSS(all_peaks, title="Distribution of peaks relative to TSS")
dev.off()

# Extract genes for GO analysis
# Human
genes_human_conserve <- unique(peakAnnoDF_human_conserve$geneId[!is.na(peakAnnoDF_human_conserve$geneId)])
genes_human_liver <- unique(peakAnnoDF_human_liver$geneId[!is.na(peakAnnoDF_human_liver$geneId)])
genes_human_pancreas <- unique(peakAnnoDF_human_pancreas$geneId[!is.na(peakAnnoDF_human_pancreas$geneId)])

# Mouse
genes_mouse_conserve <- unique(peakAnnoDF_mouse_conserve$geneId[!is.na(peakAnnoDF_mouse_conserve$geneId)])
genes_mouse_liver <- unique(peakAnnoDF_mouse_liver$geneId[!is.na(peakAnnoDF_mouse_liver$geneId)])
genes_mouse_pancreas <- unique(peakAnnoDF_mouse_pancreas$geneId[!is.na(peakAnnoDF_mouse_pancreas$geneId)])

# Create list for comparative analysis
gene_lists_human <- list(
  Conserved = genes_human_conserve,
  Liver = genes_human_liver,
  Pancreas = genes_human_pancreas
)

gene_lists_mouse <- list(
  Conserved = genes_mouse_conserve,
  Liver = genes_mouse_liver,
  Pancreas = genes_mouse_pancreas
)

# Comparative GO enrichment analysis for human
comp_go_human <- compareCluster(
  gene_lists_human,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Comparative GO enrichment analysis for mouse
comp_go_mouse <- compareCluster(
  gene_lists_mouse,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Save results to CSV
write.csv(as.data.frame(comp_go_human), "cross-tissue_results/human_GO_enrichment_results.csv")
write.csv(as.data.frame(comp_go_mouse), "cross-tissue_results/mouse_GO_enrichment_results.csv")

# Visualize comparative results
pdf("cross-tissue_results/human_GO_dotplot.pdf", width=12, height=15)
dotplot(comp_go_human, showCategory=15, title="Human Biological Processes Comparison")
dev.off()

pdf("cross-tissue_results/mouse_GO_dotplot.pdf", width=12, height=15)
dotplot(comp_go_mouse, showCategory=15, title="Mouse Biological Processes Comparison")
dev.off()


# KEGG pathway analysis
# Convert gene IDs for KEGG (might need additional steps depending on ID format)
# Human
kegg_human <- compareCluster(
  gene_lists_human,
  fun = "enrichKEGG",
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# Mouse
kegg_mouse <- compareCluster(
  gene_lists_mouse,
  fun = "enrichKEGG",
  organism = "mmu",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# Save KEGG results
write.csv(as.data.frame(kegg_human), "cross-tissue_results/human_KEGG_pathway_results.csv")
write.csv(as.data.frame(kegg_mouse), "cross-tissue_results/mouse_KEGG_pathway_results.csv")

# KEGG pathway visualization
pdf("cross-tissue_results/human_KEGG_dotplot.pdf", width=12, height=15)
dotplot(kegg_human, showCategory=15, title="Human KEGG Pathway Comparison")
dev.off()

pdf("cross-tissue_results/mouse_KEGG_dotplot.pdf", width=12, height=15)
dotplot(kegg_mouse, showCategory=15, title="Mouse KEGG Pathway Comparison")
dev.off()
