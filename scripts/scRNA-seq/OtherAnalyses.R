library(tidyverse)
library(tidyr)
library(Seurat)
library(patchwork)


harmonyEpiCells <- readRDS(file = "Harmony_integrated_epidermal_w_Merkel.RData")


#########Cell Cycle Analysis
cellCycleGenes <- read.csv("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/scRNA-seq/CellCycleGenes_Mus_musculus.csv")
s_phase <- cellCycleGenes[cellCycleGenes$phase == 'S', "geneName", drop = TRUE]
g2m_phase <- cellCycleGenes[cellCycleGenes$phase == 'G2/M', "geneName", drop = TRUE]

cellCycleGenes <- read.csv("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/scRNA-seq/Seurat_CellCycleGenes.csv", header = TRUE, row.names = 1)
s_phase <- cellCycleGenes[cellCycleGenes$Phase == 'S', "Mouse", drop = TRUE]
g2m_phase <- cellCycleGenes[cellCycleGenes$Phase == 'G2/M', "Mouse", drop = TRUE]

harmonyEpiCells <- CellCycleScoring(harmonyEpiCells, s.features = s_phase, g2m.features = g2m_phase)
DimPlot(harmonyEpiCells, group.by = 'Phase', split.by = 'genotype', pt.size = 0.6)

cellCyclePhaseInfo <- harmonyEpiCells@meta.data %>% group_by(genotype, Identity, Phase) %>%
  summarise(count = n()) %>% mutate(percent.phase = count / sum(count)) %>% ungroup()

ggplot(cellCyclePhaseInfo, aes(x = genotype, y = percent.phase, fill = Phase)) +
  geom_bar(stat = 'identity', position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c('gold', 'yellowgreen', 'darkgreen')) +
  facet_wrap(~ Identity, nrow = 1, strip.position = 'bottom') +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = '%Cells', fill = 'Phase', x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45), 
        strip.background = element_blank(), strip.placement = 'outside', strip.text.x = element_text(size = 12, angle = 60))



########Ploting heatmap for differentially regulated genes
#Ploting with DoHeatmap from Seurat package
Idents(harmonyEpiCells) <- 'seurat_clusters'
AggrExp_Seurat <- AggregateExpression(harmonyEpiCells, return.seurat = TRUE, group.by = c('genotype', 'ident'))
features_to_remove <- grepl('ENSM', scRNA_all_DE) | grepl('Rik', scRNA_all_DE) | grepl('^Gm', scRNA_all_DE)
AggrExp_Seurat$genotype_cluster <- Idents(AggrExp_Seurat)
Idents(AggrExp_Seurat) <- 'orig.ident'
AggrExp_Seurat <- RenameIdents(AggrExp_Seurat, Cluster_idents)
AggrExp_Seurat$Identity <- Idents(AggrExp_Seurat)

DoHeatmap(AggrExp_Seurat, features = scRNA_all_DE[!features_to_remove], group.by = 'genotype', lines.width = 1)

#Ploting with ComplexHeatmap package
library(ComplexHeatmap)
library(circlize)

#Get gene list for heatmap plotting
##Union
scRNA_markers_c49 <- Reduce(union, lapply(topMarkers[c(5, 10)], row.names))
scRNA_DE <- Reduce(union, Cluster_DE_geneNames[-c(5,9,10)])
features_to_plot <- union(scRNA_markers_c49, scRNA_DE)
#Intersection
scRNA_markers_c49 <- Reduce(intersect, lapply(ClusterMarkers[c(5, 10)], row.names))
scRNA_DE <- Reduce(intersect, Cluster_DE_geneNames[-c(5,9,10)])
features_to_plot <- union(scRNA_markers_c49, scRNA_DE)
#Using cluster 4 and cluster 9 markers obtained from comparing to basal or differentiated counterpart
features_to_plot <- union(scRNA_DE, Marker_C4_C9)
features_to_plot <- pseudoBulk_all_DE

features_to_remove <- grepl('ENSM', features_to_plot) | grepl('Rik', features_to_plot) | grepl('^Gm', features_to_plot)

mat_DE <- AggrExp_Seurat@assays$RNA$scale.data[features_to_plot, ]
mat_DE <- AggrExp_Seurat@assays$RNA$scale.data[VariableFeatures(harmonyEpiCells), ]
colnames(mat_DE) <- paste0(AggrExp_Seurat$genotype, '_', AggrExp_Seurat$Identity)
mat_DE <- mat_DE[, -c(5, 10, 19)]

#color_heatmap <- colorRamp2(c(-3, 0, 3), c('blue', '#EEEEEE', 'red'))
Heatmap(mat_DE, name = 'TKO vs Ctrl DE genes', col = color_heatmap, 
        cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE, column_names_side = 'top', column_names_rot = 90, row_names_side = 'left',
        column_names_gp = gpar(col = c(rep('blue', 8), rep('red', 9))),
        row_title = '1875 differentially expressed genes',
        heatmap_legend_param = list(title = 'Norm. Exp.'),
        use_raster = TRUE, raster_device = 'png'
        )


##########GSEA Analysis
#Preparing Data
library(openxlsx)

clusters <- getSheetNames("ExcelFiles/DE_genes_v3.xlsx")
for (name in clusters) {
  data <- readWorkbook("ExcelFiles/DE_genes_v3.xlsx", sheet = name)
  data <- data[,c(1,3)]
  filename <- paste0("GSEA/DE_cluster_", name, "_ranked.rnk")
  write.table(data, file = filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

data <- readWorkbook("ExcelFiles/Conserved_markers_v3.xlsx", sheet = '4')
data <- data[,c(1,3)]
write.table(data, file = "GSEA/Markers_cluster_4_ranked.rnk", quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)


Idents(harmonyEpiCells) <- "cluster_genotype"
for (name in clusters) {
  Ident_1 <- paste(name, 'TKO', sep = '_')
  Ident_2 <- paste(name, 'Ctrl', sep = '_')
  if (all(c(Ident_1, Ident_2) %in% harmonyEpiCells$cluster_genotype)) {
    cat('Finding differentially expressed genes for cluster:', name, '\n')
    DE <- FindMarkers(harmonyEpiCells, ident.1 = Ident_1, ident.2 = Ident_2, logfc.threshold = 0.3, min.pct = 0.25)
    log2FC <- DE[,2, drop = FALSE]
    filename <- paste0("GSEA/DE_all_cluster_", name, "_ranked.rnk")
    write.table(log2FC, file = filename, quote = FALSE, sep = '\t', row.names = TRUE, col.names = FALSE)
  }
}
Idents(harmonyEpiCells) <- "seurat_clusters"
rm(DE)
rm(log2FC)

marker <- FindConservedMarkers(harmonyEpiCells, ident.1 = '4', grouping.var = 'genotype', min.cells.group = 10)
log2FC <- marker[,2, drop = FALSE]
write.table(log2FC, file = "GSEA/Markers_all_cluster_4_ranked.rnk", quote = FALSE, sep = '\t', row.names = TRUE, col.names = FALSE)

marker <- FindConservedMarkers(harmonyEpiCells, ident.1 = '9', grouping.var = 'genotype', min.cells.group = 10)
log2FC <- marker[,2, drop = FALSE]
write.table(log2FC, file = "GSEA/Markers_all_cluster_9_ranked.rnk", quote = FALSE, sep = '\t', row.names = TRUE, col.names = FALSE)


#######Test if methylation predicts gene expression change
methyl_flank <- read.xlsx("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/Methyl_annotation_allRegions_allGenes.xlsx",
                                      sheet = "flank10k")
methyl_geneBody <- read.xlsx("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/Methyl_annotation_allRegions_allGenes.xlsx",
                                      sheet = "geneBody")
methyl_promoter <- read.xlsx("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/Methyl_annotation_allRegions_allGenes.xlsx",
                          sheet = "promoter")

methyl_associated_genes <- subset(methyl_promoter, NumPeaks > 0)
methyl_associated_genes <- subset(methyl_geneBody, MaxPeak > 3)
methyl_associated_genes <- subset(methyl_geneBody, MeanIntensity > 0.5)
methyl_associated_genes <- subset(methyl_flank, MaxPeak > 3)
methyl_associated_genes <- subset(methyl_flank, MeanIntensity > 0.5)

mat_DE <- AggrExp_Seurat@assays$RNA$scale.data
colnames(mat_DE) <- paste0(AggrExp_Seurat$genotype, '_', AggrExp_Seurat$Identity)
keep <- rownames(mat_DE) %in% methyl_associated_genes$ensembl | rownames(mat_DE) %in% methyl_associated_genes$name
mat_DE <- mat_DE[keep, -c(5, 10, 19)]
#mat_DE <- mat_DE[rownames(mat_DE) != 'Cd300ld', ]

Heatmap(mat_DE, name = 'TKO vs Ctrl DE genes', col = color_heatmap, 
        cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE, cluster_columns = FALSE,
        column_names_side = 'top', column_names_rot = 90, row_names_side = 'left',
        column_names_gp = gpar(col = c(rep('blue', 8), rep('red', 9))),
        heatmap_legend_param = list(title = 'Norm. Exp.'),
        use_raster = TRUE, raster_device = 'png'
)


DEgenes_Up <- lapply(Cluster_DE_genes, function(df) rownames(df[df$avg_log2FC > 0,]))
Up_genes_combined <- Reduce(union, DEgenes_Up[-c(5,9,10)])
ggvenn(list('Predicted_up' = rownames(mat_DE), 'scRNAseq_up' = Up_genes_combined), 
       stroke_size = 1, set_name_size = 6, text_size = 5) + expand_limits(y = 2)

rm(list = ls()[grep('^methyl', ls())])





