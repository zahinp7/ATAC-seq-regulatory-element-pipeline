#############################################################
# Functional annotation of conserved and specific regions   #
# between human and mouse in liver and pancreas            #
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

# Set file paths (adjust as needed)
args <- commandArgs(trailingOnly = TRUE)
input_file1 <- args[1]
input_file2 <- args[2]
input_file3 <- args[3]
input_file4 <- args[4]
input_file5 <- args[5]
input_file6 <- args[6]

# Import BED files
# Conservation and specificity files
liver_conserved <- readPeakFile(input_file1)
human_liver_specific <- readPeakFile(input_file2)
mouse_liver_specific <- readPeakFile(input_file3)

pancreas_conserved <- readPeakFile(input_file4)
human_pancreas_specific <- readPeakFile(input_file5)
mouse_pancreas_specific <- readPeakFile(input_file6)

# Set up annotation databases
txdb_human <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb_mouse <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Annotate peaks
# Use human TxDb for human regions and conserved regions (represented in human coordinates)
liver_conserved_anno <- annotatePeak(liver_conserved, tssRegion=c(-3000, 3000),
                                     TxDb=txdb_human, annoDb="org.Hs.eg.db")
human_liver_specific_anno <- annotatePeak(human_liver_specific, tssRegion=c(-3000, 3000),
                                          TxDb=txdb_human, annoDb="org.Hs.eg.db")
pancreas_conserved_anno <- annotatePeak(pancreas_conserved, tssRegion=c(-3000, 3000),
                                        TxDb=txdb_human, annoDb="org.Hs.eg.db")
human_pancreas_specific_anno <- annotatePeak(human_pancreas_specific, tssRegion=c(-3000, 3000),
                                             TxDb=txdb_human, annoDb="org.Hs.eg.db")

# Use mouse TxDb for mouse-specific regions
mouse_liver_specific_anno <- annotatePeak(mouse_liver_specific, tssRegion=c(-3000, 3000),
                                          TxDb=txdb_mouse, annoDb="org.Mm.eg.db")
mouse_pancreas_specific_anno <- annotatePeak(mouse_pancreas_specific, tssRegion=c(-3000, 3000),
                                             TxDb=txdb_mouse, annoDb="org.Mm.eg.db")

# Create a directory for results
dir.create("cross-species_results", showWarnings = FALSE)

# Create annotation lists for visualization
liver_anno_list <- list(
  "Conserved" = liver_conserved_anno,
  "Human_Specific" = human_liver_specific_anno,
  "Mouse_Specific" = mouse_liver_specific_anno
)

pancreas_anno_list <- list(
  "Conserved" = pancreas_conserved_anno,
  "Human_Specific" = human_pancreas_specific_anno,
  "Mouse_Specific" = mouse_pancreas_specific_anno
)

# Visualize genomic feature distribution
# Liver genomic annotations
pdf("cross-species_results/liver_genomic_features_bar.pdf", width=10, height=6)
plotAnnoBar(liver_anno_list)
dev.off()

# Pancreas genomic annotations
pdf("cross-species_results/pancreas_genomic_features_bar.pdf", width=10, height=6)
plotAnnoBar(pancreas_anno_list)
dev.off()

# Convert to data frames
liver_conserved_df <- as.data.frame(liver_conserved_anno)
human_liver_specific_df <- as.data.frame(human_liver_specific_anno)
mouse_liver_specific_df <- as.data.frame(mouse_liver_specific_anno)
pancreas_conserved_df <- as.data.frame(pancreas_conserved_anno)
human_pancreas_specific_df <- as.data.frame(human_pancreas_specific_anno)
mouse_pancreas_specific_df <- as.data.frame(mouse_pancreas_specific_anno)

# Extract genes for functional analysis
liver_conserved_genes <- unique(liver_conserved_df$geneId[!is.na(liver_conserved_df$geneId)])
human_liver_specific_genes <- unique(human_liver_specific_df$geneId[!is.na(human_liver_specific_df$geneId)])
mouse_liver_specific_genes <- unique(mouse_liver_specific_df$geneId[!is.na(mouse_liver_specific_df$geneId)])
pancreas_conserved_genes <- unique(pancreas_conserved_df$geneId[!is.na(pancreas_conserved_df$geneId)])
human_pancreas_specific_genes <- unique(human_pancreas_specific_df$geneId[!is.na(human_pancreas_specific_df$geneId)])
mouse_pancreas_specific_genes <- unique(mouse_pancreas_specific_df$geneId[!is.na(mouse_pancreas_specific_df$geneId)])

# Create gene lists with consistent structure across all categories
liver_gene_list <- list(
  "Conserved" = liver_conserved_genes,
  "Human_Specific" = human_liver_specific_genes
)

pancreas_gene_list <- list(
  "Conserved" = pancreas_conserved_genes,
  "Human_Specific" = human_pancreas_specific_genes
)

# Mouse-specific gene lists
mouse_liver_gene_list <- list(
  "Mouse_Specific" = mouse_liver_specific_genes
)

mouse_pancreas_gene_list <- list(
  "Mouse_Specific" = mouse_pancreas_specific_genes
)

###################### GO ENRICHMENT ANALYSIS ######################
# Human genes (conserved and human-specific)
liver_go_human <- compareCluster(
  liver_gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

pancreas_go_human <- compareCluster(
  pancreas_gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Mouse-specific genes
liver_go_mouse <- compareCluster(
  mouse_liver_gene_list,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

pancreas_go_mouse <- compareCluster(
  mouse_pancreas_gene_list,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Save GO results
write.csv(as.data.frame(liver_go_human), "cross-species_results/liver_GO_human_results.csv")
write.csv(as.data.frame(liver_go_mouse), "cross-species_results/liver_GO_mouse_results.csv")
write.csv(as.data.frame(pancreas_go_human), "cross-species_results/pancreas_GO_human_results.csv")
write.csv(as.data.frame(pancreas_go_mouse), "cross-species_results/pancreas_GO_mouse_results.csv")

###################### COMBINED GO VISUALIZATION ######################
# Extract data from GO results for combined visualization
extract_top_GO <- function(go_result, n=10) {
  if (is.null(go_result) || nrow(as.data.frame(go_result)) == 0) {
    return(data.frame())
  }
  
  df <- as.data.frame(go_result)
  
  # Get top n categories for each cluster
  top_df <- data.frame()
  for (cluster in unique(df$Cluster)) {
    cluster_df <- df[df$Cluster == cluster, ]
    cluster_df <- cluster_df[order(cluster_df$p.adjust), ]
    cluster_df <- cluster_df[1:min(n, nrow(cluster_df)), ]
    top_df <- rbind(top_df, cluster_df)
  }
  
  return(top_df)
}

# Extract top GO data
liver_go_human_data <- extract_top_GO(liver_go_human)
liver_go_mouse_data <- extract_top_GO(liver_go_mouse)
pancreas_go_human_data <- extract_top_GO(pancreas_go_human)
pancreas_go_mouse_data <- extract_top_GO(pancreas_go_mouse)

# Add source information
if (nrow(liver_go_human_data) > 0) liver_go_human_data$Source <- "Human"
if (nrow(liver_go_mouse_data) > 0) liver_go_mouse_data$Source <- "Mouse" 
if (nrow(pancreas_go_human_data) > 0) pancreas_go_human_data$Source <- "Human"
if (nrow(pancreas_go_mouse_data) > 0) pancreas_go_mouse_data$Source <- "Mouse"

# Combine data for visualization
liver_go_combined <- rbind(liver_go_human_data, liver_go_mouse_data)
pancreas_go_combined <- rbind(pancreas_go_human_data, pancreas_go_mouse_data)

# Create custom GO dotplots for liver
if (nrow(liver_go_combined) > 0) {
  # Truncate long descriptions
  liver_go_combined$Description <- substr(liver_go_combined$Description, 1, 50)
  
  # Create combined plot for liver
  pdf("cross-species_results/liver_GO_combined_dotplot.pdf", width=14, height=12)
  p <- ggplot(liver_go_combined, 
              aes(x=Count, y=Description, size=Count, color=-log10(p.adjust))) +
    geom_point() +
    facet_grid(Cluster ~ ., scales="free_y", space="free") +
    scale_color_continuous(name="-log10(p.adjust)") +
    scale_size_continuous(name="Gene Count") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size=8),
      strip.text = element_text(size=10, face="bold"),
      plot.title = element_text(hjust=0.5, size=12, face="bold")
    ) +
    labs(title="Liver - Biological Processes (Combined)", x="Gene Count", y="")
  print(p)
  dev.off()
}

# Create custom GO dotplots for pancreas
if (nrow(pancreas_go_combined) > 0) {
  # Truncate long descriptions
  pancreas_go_combined$Description <- substr(pancreas_go_combined$Description, 1, 50)
  
  # Create combined plot for pancreas
  pdf("cross-species_results/pancreas_GO_combined_dotplot.pdf", width=14, height=12)
  p <- ggplot(pancreas_go_combined, 
              aes(x=Count, y=Description, size=Count, color=-log10(p.adjust))) +
    geom_point() +
    facet_grid(Cluster ~ ., scales="free_y", space="free") +
    scale_color_continuous(name="-log10(p.adjust)") +
    scale_size_continuous(name="Gene Count") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size=8),
      strip.text = element_text(size=10, face="bold"),
      plot.title = element_text(hjust=0.5, size=12, face="bold")
    ) +
    labs(title="Pancreas - Biological Processes (Combined)", x="Gene Count", y="")
  print(p)
  dev.off()
}

###################### KEGG PATHWAY ANALYSIS ######################
# Human genes (conserved and human-specific)
liver_kegg_human <- compareCluster(
  liver_gene_list,
  fun = "enrichKEGG",
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

pancreas_kegg_human <- compareCluster(
  pancreas_gene_list,
  fun = "enrichKEGG",
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# Mouse-specific genes
liver_kegg_mouse <- compareCluster(
  mouse_liver_gene_list,
  fun = "enrichKEGG",
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

pancreas_kegg_mouse <- compareCluster(
  mouse_pancreas_gene_list,
  fun = "enrichKEGG",
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# Save KEGG results
write.csv(as.data.frame(liver_kegg_human), "cross-species_results/liver_KEGG_human_results.csv")
write.csv(as.data.frame(liver_kegg_mouse), "cross-species_results/liver_KEGG_mouse_results.csv")
write.csv(as.data.frame(pancreas_kegg_human), "cross-species_results/pancreas_KEGG_human_results.csv")
write.csv(as.data.frame(pancreas_kegg_mouse), "cross-species_results/pancreas_KEGG_mouse_results.csv")

###################### COMBINED KEGG VISUALIZATION ######################
# Extract data from KEGG results
extract_top_KEGG <- function(kegg_result, n=10) {
  if (is.null(kegg_result) || nrow(as.data.frame(kegg_result)) == 0) {
    return(data.frame())
  }
  
  df <- as.data.frame(kegg_result)
  
  # Get top n categories for each cluster
  top_df <- data.frame()
  for (cluster in unique(df$Cluster)) {
    cluster_df <- df[df$Cluster == cluster, ]
    cluster_df <- cluster_df[order(cluster_df$p.adjust), ]
    cluster_df <- cluster_df[1:min(n, nrow(cluster_df)), ]
    top_df <- rbind(top_df, cluster_df)
  }
  
  return(top_df)
}

# Extract top KEGG data
liver_kegg_human_data <- extract_top_KEGG(liver_kegg_human)
liver_kegg_mouse_data <- extract_top_KEGG(liver_kegg_mouse)
pancreas_kegg_human_data <- extract_top_KEGG(pancreas_kegg_human)
pancreas_kegg_mouse_data <- extract_top_KEGG(pancreas_kegg_mouse)

# Add source information
if (nrow(liver_kegg_human_data) > 0) liver_kegg_human_data$Source <- "Human"
if (nrow(liver_kegg_mouse_data) > 0) liver_kegg_mouse_data$Source <- "Mouse"
if (nrow(pancreas_kegg_human_data) > 0) pancreas_kegg_human_data$Source <- "Human"
if (nrow(pancreas_kegg_mouse_data) > 0) pancreas_kegg_mouse_data$Source <- "Mouse"

# Combine data for visualization
liver_kegg_combined <- rbind(liver_kegg_human_data, liver_kegg_mouse_data)
pancreas_kegg_combined <- rbind(pancreas_kegg_human_data, pancreas_kegg_mouse_data)

# Create custom KEGG dotplots for liver
if (nrow(liver_kegg_combined) > 0) {
  # Truncate long descriptions
  liver_kegg_combined$Description <- substr(liver_kegg_combined$Description, 1, 50)
  
  # Create combined plot for liver
  pdf("cross-species_results/liver_KEGG_combined_dotplot.pdf", width=14, height=12)
  p <- ggplot(liver_kegg_combined, 
              aes(x=Count, y=Description, size=Count, color=-log10(p.adjust))) +
    geom_point() +
    facet_grid(Cluster ~ ., scales="free_y", space="free") +
    scale_color_continuous(name="-log10(p.adjust)") +
    scale_size_continuous(name="Gene Count") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size=8),
      strip.text = element_text(size=10, face="bold"),
      plot.title = element_text(hjust=0.5, size=12, face="bold")
    ) +
    labs(title="Liver - KEGG Pathways (Combined)", x="Gene Count", y="")
  print(p)
  dev.off()
}

# Create custom KEGG dotplots for pancreas
if (nrow(pancreas_kegg_combined) > 0) {
  # Truncate long descriptions
  pancreas_kegg_combined$Description <- substr(pancreas_kegg_combined$Description, 1, 50)
  
  # Create combined plot for pancreas
  pdf("cross-species_results/pancreas_KEGG_combined_dotplot.pdf", width=14, height=12)
  p <- ggplot(pancreas_kegg_combined, 
              aes(x=Count, y=Description, size=Count, color=-log10(p.adjust))) +
    geom_point() +
    facet_grid(Cluster ~ ., scales="free_y", space="free") +
    scale_color_continuous(name="-log10(p.adjust)") +
    scale_size_continuous(name="Gene Count") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size=8),
      strip.text = element_text(size=10, face="bold"),
      plot.title = element_text(hjust=0.5, size=12, face="bold")
    ) +
    labs(title="Pancreas - KEGG Pathways (Combined)", x="Gene Count", y="")
  print(p)
  dev.off()
}
