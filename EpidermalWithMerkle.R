library(tidyverse)
library(tidyr)
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)

filteredCells <- readRDS(file = "Filtered cells.RData")
EpidermalCells <- subset(filteredCells, idents = c('0', '1', '2', '3', '4', '7', '8', 'Merkel cells'), invert = FALSE)

#Renormalize and cluster
EpidermalCells <- NormalizeData(EpidermalCells, normalization.method = "LogNormalize", scale.factor = 10000)
EpidermalCells <- FindVariableFeatures(EpidermalCells, selection.method = "vst", nfeatures = 2000)
all.features <- rownames(EpidermalCells)
EpidermalCells <- ScaleData(EpidermalCells, features = all.features)

EpidermalCells <- RunPCA(EpidermalCells, features = VariableFeatures(EpidermalCells), reduction.name = 'pca')
ElbowPlot(EpidermalCells)
EpidermalCells <- FindNeighbors(EpidermalCells, dims = 1:20)
EpidermalCells <- FindClusters(EpidermalCells, resolution = 0.25)
EpidermalCells <- RunUMAP(EpidermalCells, dims = 1:20, reduction.name = 'umap_no_int')
DimPlot(EpidermalCells, reduction = "umap_no_int", split.by = "genotype", pt.size = 0.6, label = TRUE)

#Integration analyis with harmony
library(harmony)

harmonyEpiCells <- IntegrateLayers(EpidermalCells, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony")
harmonyEpiCells[['RNA']] <- JoinLayers(harmonyEpiCells[['RNA']])

harmonyEpiCells <- FindNeighbors(harmonyEpiCells, reduction = 'harmony', dims = 1:20)
harmonyEpiCells <- FindClusters(harmonyEpiCells, resolution = 0.2)
harmonyEpiCells <- RunUMAP(harmonyEpiCells, reduction = 'harmony', dims = 1:20, reduction.name = 'umap')

DimPlot(harmonyEpiCells, reduction = 'umap', split.by = 'genotype', pt.size = 0.6, label.size = 5, label = TRUE)

# Plot the percentage of each cluster
cluster_count <- data.frame(sample = harmonyEpiCells$genotype, cluster = Idents(harmonyEpiCells)) %>% 
  group_by(sample, cluster) %>% summarise(count = n(), .groups = 'drop')
additional <- data.frame(sample = c('Ctrl', 'Ctrl', 'TKO'),
                         cluster = c('AberrBasal', 'AberrDiff', 'Merkel'),
                         count = c(0, 0, 0))
cluster_count <- rbind(cluster_count, additional)
sample_sum <- cluster_count %>% group_by(sample) %>% summarise(sum_per_sample = sum(count), .groups = 'drop')
cluster_pct <- cluster_count %>% left_join(sample_sum, by = 'sample') %>%
  mutate(percentage = count/sum_per_sample) %>%
  mutate(sample = factor(sample, levels = c('TKO', 'Ctrl')))
cluster_pct$cluster <- fct_rev(cluster_pct$cluster)

cluster_colors <- rev(RColorBrewer::brewer.pal(10, 'Set3'))
cluster_colors[2] <- rgb(217, 187, 140, maxColorValue = 255)
cluster_colors[3] <- rgb(252, 214, 17, maxColorValue = 255)

ggplot(cluster_pct, aes(x = sample, y = percentage, fill = cluster)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = cluster_colors) +
  #scale_fill_brewer(palette = 'Set3', direction = -1) +
  scale_y_continuous(expand = c(0.01, 0)) +
  coord_flip() +
  labs(x = NULL, y = 'Percentage', fill = 'Clusters')

cluster_pct$sample <- factor(cluster_pct$sample, levels = c('Ctrl', 'TKO'))
cluster_levels <- c('Basal_1', 'Basal_2', 'Basal_3', 'AberrBasal', 'Diff_1', 'Diff_2', 'Diff_3', 'AberrDiff', 'HairFollicle', 'Merkel')
cluster_pct$cluster <- factor(cluster_pct$cluster, levels = cluster_levels)
ggplot(cluster_pct, aes(x = sample, y = percentage, fill = cluster)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = rev(cluster_colors)) +
  #scale_fill_brewer(palette = 'Set3', direction = -1) +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(x = NULL, y = 'Percentage', fill = 'Clusters')

ggplot(cluster_pct, aes(x = cluster, y = percentage, fill = sample)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = sprintf("%.2f", percentage)), position = position_dodge(width = 0.9), vjust = -0.5, size = 2.5) +
  scale_fill_manual(values = c('darkgreen', 'gold')) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 0.45)) +
  labs(x = NULL, y = 'Percentage', fill = 'Genotype') +
  theme(axis.title.y = element_text(size = 12))

rm(cluster_count, sample_sum, cluster_pct)

# Check the expression of a few genes
VlnPlot(harmonyEpiCells, features = c('Krt14', 'Itga6', 'Krt10', 'Dsg1a', 'Lor', 'Mki67', 'Top2a', 'Sox9', 'Nfatc1'), ncol = 3, alpha = 0)
FeaturePlot(harmonyEpiCells, features = 'Krt14', split.by = 'genotype', pt.size = 0.6)
VlnPlot(harmonyEpiCells, features = c('Mki67', 'Top2a'), split.by = 'genotype', alpha = 0.3, cols = c('lightblue', 'coral'), ncol = 2) +
  theme(legend.position = 'right')
VlnPlot(harmonyEpiCells, features = c('Col28a1', 'Hopx', 'Cd28', 'Gria3'), split.by = 'genotype', alpha = 0, cols = c('lightblue', 'coral'), ncol = 2) +
  theme(legend.position = 'right')

#Assign cluster identities
Cluster_idents <- c('0' = 'Basal_Qscnt', '1' = 'Diff_main', '2' = 'Basal_Proli', '3' = 'Hair_follicle',
                    '4' = 'TKO_new', '5' = 'Hair_placode', '6' = 'Diff_sub1', '7' = 'Diff_sub2',
                    '8' = 'Merkel', '9' = 'TKO_new_2')
harmonyEpiCells <- RenameIdents(harmonyEpiCells, Cluster_idents)
harmonyEpiCells$Identity <- Idents(harmonyEpiCells)
Idents(harmonyEpiCells) <- "seurat_clusters"
Idents(harmonyEpiCells) <- "Identity"
DimPlot(harmonyEpiCells, reduction = 'umap', split.by = 'genotype', group.by = 'Identity', pt.size = 0.6, label.size = 5, label = FALSE)


#Find markers for each cluster
library(multtest)
library(metap)
library(openxlsx)

ClusterMarkers <- list()
clusters <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
for (name in c('0', '1', '2', '3', '5', '6', '7')) {
  cat('Finding markers for cluster:', name, '\n')
  markers <- FindConservedMarkers(harmonyEpiCells, ident.1 = name, grouping.var = 'genotype', min.cells.group = 10)
  keep1 <- markers$TKO_p_val_adj < 0.05 & abs(markers$TKO_avg_log2FC) > 1
  keep2 <- markers$Ctrl_p_val_adj < 0.05 & abs(markers$Ctrl_avg_log2FC) > 1
  ClusterMarkers[[name]] <- subset(markers, keep1 | keep2)
}

for (name in c('4', '8', '9')) {
  cat('Finding markers for cluster:', name, '\n')
  markers <- FindConservedMarkers(harmonyEpiCells, ident.1 = name, grouping.var = 'genotype', min.cells.group = 10)
  if (is.null(markers$TKO_avg_log2FC)) keep <- markers$Ctrl_p_val_adj < 0.05 & abs(markers$Ctrl_avg_log2FC) > 1
  else keep <- markers$TKO_p_val_adj < 0.05 & abs(markers$TKO_avg_log2FC) > 1
  
  ClusterMarkers[[name]] <- subset(markers, keep)
}

wb <- createWorkbook()
for (name in clusters) {
  cat('Writing data to excel for cluster:', name, '\n')
  addWorksheet(wb, name)
  writeData(wb, name, ClusterMarkers[[name]], rowNames = TRUE)
}
saveWorkbook(wb, 'ExcelFiles/Conserved_markers_v3.xlsx', overwrite = TRUE)
rm(wb)
rm(markers)

topMarkers <- lapply(ClusterMarkers, function(marker){
  keep_1 <- FALSE
  keep_2 <- FALSE
  if ('TKO_p_val_adj' %in% names(marker)) {
    keep_1 <- (marker$TKO_pct.1 - marker$TKO_pct.2 > 0.4) | marker$TKO_avg_log2FC > 3
  }
  if ('Ctrl_p_val_adj' %in% names(marker)) {
    keep_2 <- marker$Ctrl_avg_log2FC > 3 | (marker$Ctrl_pct.1 - marker$Ctrl_pct.2 > 0.4)
  }
  marker <- subset(marker, subset = keep_1 | keep_2)
})

##########################################################
#Differential analysis bewteen TKO and Ctrl
harmonyEpiCells$geno_cluster <- paste(harmonyEpiCells$genotype, harmonyEpiCells$seurat_clusters, sep = '_')

Idents(harmonyEpiCells) <- "geno_cluster"
Cluster_DE_genes <- list()
clusters <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
for (name in clusters) {
  Ident_1 <- paste('TKO', name, sep = '_')
  Ident_2 <- paste('Ctrl', name, sep = '_')
  cat('Finding differentially expressed genes for cluster:', name, '\n')
  DE <- FindMarkers(harmonyEpiCells, ident.1 = Ident_1, ident.2 = Ident_2, logfc.threshold = 0.3, min.pct = 0.25)
  Cluster_DE_genes[[name]] <- subset(DE, abs(avg_log2FC) > 1 & p_val_adj < 0.05)
}

DE <- FindMarkers(harmonyEpiCells, ident.1 = 'TKO_4', ident.2 = c('Ctrl_0', 'Ctrl_2'), logfc.threshold = 0.3, min.pct = 0.25)
Cluster_DE_genes[['4']] <- subset(DE, abs(avg_log2FC) > 1 & p_val_adj < 0.05)

DE <- FindMarkers(harmonyEpiCells, ident.1 = 'TKO_9', ident.2 = c('Ctrl_1', 'Ctrl_6', 'Ctrl_7'), logfc.threshold = 0.3, min.pct = 0.25)
Cluster_DE_genes[['9']] <- subset(DE, abs(avg_log2FC) > 1 & p_val_adj < 0.05)

Idents(harmonyEpiCells) <- "seurat_clusters"

wb <- createWorkbook()
for (name in clusters) {
  cat('Writing data to excel for cluster:', name, '\n')
  addWorksheet(wb, name)
  writeData(wb, name, Cluster_DE_genes[[name]], rowNames = TRUE)
}
saveWorkbook(wb, 'ExcelFiles/DE_genes_v4.xlsx', overwrite = TRUE)
addGeneInfo('ExcelFiles/DE_genes_v4.xlsx', ensembl = Ens_mm39_110, ID_type = 'SYMBOL')
rm(wb)
rm(DE)

############################################################
#Co-analyzing with RNA-seq DE genes
library(ggvenn)

RNAseq_DE <- read.csv(file = "TKOvsCtrl.csv", header = TRUE, row.names = 1)
RNAseq_DE <- subset(RNAseq_DE, abs(RNAseq_DE$log2FoldChange) > 1 & RNAseq_DE$pvalue < 0.01)
RNAseq_DE <- na.omit(RNAseq_DE)

RNAseq_DE$Ensembl <- row.names(RNAseq_DE)
row.names(RNAseq_DE) <- RNAseq_DE$symbol

Cluster_DE_geneNames <- lapply(Cluster_DE_genes, row.names)
scRNA_all_DE <- Reduce(union, Cluster_DE_geneNames[-c(5,9,10)])

####Subset to remove clusters with < 10 cells for plots
Idents(harmonyEpiCells) <- 'geno_cluster'
harmonyCells <- subset(harmonyEpiCells, idents = c('Ctrl_4', 'Ctrl_9', 'TKO_8'), invert = TRUE)
Idents(harmonyCells) <- 'Identity'
rm(harmonyEpiCells)

features_to_plot <- c('Krt14', 'Krt10', 'Col28a1', 'Kyat3', 'Gria3', 'Cd28', 'Ly6e', 'Tcfl5', 'Col4a6', 'Lama1')
VlnPlot(harmonyCells, features = features_to_plot, split.by = 'genotype', alpha = 0, cols = c('lightblue', 'coral'), ncol = 2) +
  theme(legend.position = 'right')
#Merkel cell related features
VlnPlot(harmonyEpiCells, features = c('Krt17', 'Krt8', 'Krt18', 'Krt20'), split.by = 'genotype', alpha = 0, ncol = 2)
#Downregulated extracelluar matrix regulated genes
VlnPlot(harmonyCells, features = c('Postn', 'Stfa3'), split.by = 'genotype', alpha = 0, cols = c('lightblue', 'coral'), ncol = 1)


#Closing and Opening
saveRDS(harmonyEpiCells, file = "Harmony_integrated_epidermal_w_Merkel.RData")
rm(harmonyEpiCells)
save.image(file = "Workspace_EpidermalwMerkel.RData")

harmonyEpiCells <- readRDS(file = "Harmony_integrated_epidermal_w_Merkel.RData")

