library(Seurat)
library(tidyverse)

harmonyEpiCells <- readRDS(file = "Harmony_integrated_epidermal_w_Merkel.RData")

Sub_HairFollicle <- subset(harmonyEpiCells, ident = c('HairFollicle'), invert = FALSE)
DimPlot(Sub_HairFollicle, reduction = 'umap', split.by = 'genotype', pt.size = 3.6, label.size = 5,
        label = FALSE, raster = TRUE)

Sub_HairFollicle <- FindNeighbors(Sub_HairFollicle, reduction = 'harmony', dims = 1:20)
Sub_HairFollicle <- FindClusters(Sub_HairFollicle, resolution = 0.2)
Sub_HairFollicle <- RunUMAP(Sub_HairFollicle, reduction = 'harmony', dims = 1:20, reduction.name = 'umap')

DimPlot(Sub_HairFollicle, reduction = 'umap', split.by = 'genotype', pt.size = 0.6, label.size = 5, label = TRUE)

cluster_count <- data.frame(sample = Sub_HairFollicle$genotype, cluster = Idents(Sub_HairFollicle)) %>% 
  group_by(sample, cluster) %>% summarise(count = n(), .groups = 'drop')
sample_sum <- cluster_count %>% group_by(sample) %>% summarise(sum_per_sample = sum(count), .groups = 'drop')
cluster_pct <- cluster_count %>% left_join(sample_sum, by = 'sample') %>%
  mutate(percentage = count/sum_per_sample)
ggplot(cluster_pct, aes(x = sample, y = percentage, fill = cluster)) +
  geom_bar(stat = 'identity') +
  scale_y_continuous(expand = c(0.01, 0)) +
  #scale_fill_brewer(palette = 'RdYlBu') +
  labs(x = NULL, y = 'Percentage', fill = 'Clusters')

rm(cluster_count, sample_sum, cluster_pct)

saveRDS(Sub_HairFollicle, file = "Harmony_HairFolllicle_Subpopulation.RData")
Sub_HairFollicle <- readRDS(file = "Harmony_HairFolllicle_Subpopulation.RData")

#################

FeaturePlot(Sub_HairFollicle, features = c('Gli2', 'Sox9'), split.by = 'genotype')
FeaturePlot(Sub_HairFollicle, features = c('Dkk2', 'Lef1'), split.by = 'genotype')
FeaturePlot(Sub_HairFollicle, features = c('Lhx2', 'Runx1', 'Sox9', 'Nfatc1'))
FeaturePlot(Sub_HairFollicle, features = c('Eda', 'Edar', 'Nfkb1'), split.by = 'genotype')
FeaturePlot(Sub_HairFollicle, features = c('Shh', 'Hoxc13', 'Gli1', 'Msx1', 'Msx2', 'Lgr5', 'Lgr6', 'Pdgfra', 'Runx1'))
FeaturePlot(Sub_HairFollicle, features = c('Wls', 'Lef1', 'Wnt10b', 'Ctnnb1', 'Dkk4'))

VlnPlot(Sub_HairFollicle, features = c('Shh', 'Mki67', 'Nfatc1', 'Krt79'))

#Dotplot showing marker genes validating cell identities
marker_genes <- c('Shh', 'Edar', 'Fgf20', 'Msx2', 'Dkk4', 'Lef1', 'Mki67', 'Top2a', 'Lhx2', 'Lrig1', 'Sox9', 'Dlg2', 'Krt10')

DotPlot(Sub_HairFollicle, features = marker_genes, cols = 'RdBu') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7),
        axis.title = element_blank())

# Find markers for each cluster ########
library(multtest)
library(metap)
library(openxlsx)

HF_Markers <- list()
for (name in c('0', '1', '2', '3')) {
  cat('Finding markers for cluster:', name, '\n')
  markers <- FindConservedMarkers(Sub_HairFollicle, ident.1 = name, grouping.var = 'genotype', min.cells.group = 10)
  keep1 <- markers$TKO_p_val_adj < 0.05 & abs(markers$TKO_avg_log2FC) > 1
  keep2 <- markers$Ctrl_p_val_adj < 0.05 & abs(markers$Ctrl_avg_log2FC) > 1
  HF_Markers[[name]] <- subset(markers, keep1 | keep2)
}

wb <- createWorkbook()
for (name in c('0', '1', '2', '3')) {
  cat('Writing data to excel for cluster:', name, '\n')
  addWorksheet(wb, name)
  writeData(wb, name, HF_Markers[[name]], rowNames = TRUE)
}
saveWorkbook(wb, 'ExcelFiles/HFsub_markers_v1.xlsx', overwrite = TRUE)
rm(wb, markers)

