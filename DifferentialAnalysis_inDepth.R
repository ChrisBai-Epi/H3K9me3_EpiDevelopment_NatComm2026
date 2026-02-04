library(tidyverse)
library(tidyr)
library(ComplexHeatmap)
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(openxlsx)

harmonyEpiCells <- readRDS(file = "Harmony_integrated_epidermal_w_Merkel.RData")
harmonyEpiCells$geno_cluster <- paste0(harmonyEpiCells$genotype, '_', harmonyEpiCells$seurat_clusters)

Cluster_idents <- c('0' = 'Basal_Qscnt', '1' = 'Diff_main', '2' = 'Basal_Proli',
                    '3' = 'Hair_follicle', '4' = 'TKO_new', '5' = 'Hair_placode',
                    '6' = 'Diff_sub1', '7' = 'Diff_sub2', '8' = 'Merkel',
                    '9' = 'TKO_new_2')

Idents(harmonyEpiCells) <- 'Identity'

##Differential Analysis between clusters
#Between TKO population and Ctrl populations
Idents(harmonyEpiCells) <- 'genotype'
DE <- FindMarkers(harmonyEpiCells, ident.1 = 'TKO', ident.2 = 'Ctrl', logfc.threshold = 0.3, min.pct = 0.25)
keep <- DE$p_val_adj < 0.05 & (abs(DE$avg_log2FC) > 1 | abs(DE$pct.1 - DE$pct.2) > 0.6)
topDE <- subset(DE, keep)

mat_DE <- AggrExp_Seurat@assays$RNA$scale.data[rownames(topDE),]
mat_DE <- AggrExp_Seurat@assays$RNA$scale.data[intersect(RNAseq_E16_DEs$GeneID, all.features),]
colnames(mat_DE) <- paste0(AggrExp_Seurat$genotype, '_', AggrExp_Seurat$Identity)
mat_DE <- mat_DE[, -c(5, 10, 19)]

Heatmap(mat_DE, name = 'TKO vs Ctrl DE genes', col = color_heatmap, 
        cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE, column_names_side = 'top', column_names_rot = 90,
        column_names_gp = gpar(col = c(rep('blue', 8), rep('red', 9))),
        row_title = '225 differentially expressed genes',
        heatmap_legend_param = list(title = 'Norm. Exp.')
)

topDE_up <- topDE[topDE$avg_log2FC > 0, ]
go_DE <- enrichGO(rownames(topDE_up), universe = all.features,
                  OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go_DE, showCategory = 20, font.size = 9)

topDE_down <- topDE[topDE$avg_log2FC < 0, ]
go_DE <- enrichGO(rownames(topDE_down), universe = all.features,
                      OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go_DE, showCategory = 10, font.size = 9) + labs(title = '88 Down-regulated genes')

write.xlsx(topDE, file = 'TKO_all_vs_Ctrl_all.xlsx', rowNames = TRUE)


#Between TKO basal population and Ctrl basal population
Idents(harmonyEpiCells) <- 'geno_cluster'
DE <- FindMarkers(harmonyEpiCells, ident.1 = c('TKO_0', 'TKO_2', 'TKO_4'), ident.2 = c('Ctrl_0', 'Ctrl_2'), logfc.threshold = 0.3, min.pct = 0.25)
keep <- DE$p_val_adj < 0.05 & (DE$avg_log2FC > 1 | abs(DE$pct.1 - DE$pct.2) > 0.6 | DE$avg_log2FC < -1)
topDE <- subset(DE, keep)

mat_DE <- AggrExp_Seurat@assays$RNA$scale.data[rownames(topDE),]
colnames(mat_DE) <- paste0(AggrExp_Seurat$genotype, '_', AggrExp_Seurat$Identity)
mat_DE <- mat_DE[, -c(5, 10, 19)]

Heatmap(mat_DE, name = 'TKO vs Ctrl DE genes', col = color_heatmap, 
        cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE, column_names_side = 'top', column_names_rot = 90,
        column_names_gp = gpar(col = c(rep('blue', 8), rep('red', 9))),
        row_title = '314 differentially expressed genes',
        heatmap_legend_param = list(title = 'Norm. Exp.')
)

topDE_up <- topDE[topDE$avg_log2FC > 0, ]
go_DE_up <- enrichGO(rownames(topDE_up), universe = all.features,
                  OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go_DE_up, showCategory = 10, font.size = 9) + labs(title = '194 Up-regulated genes')

topDE_down <- topDE[topDE$avg_log2FC < 0, ]
go_DE_down <- enrichGO(rownames(topDE_down), universe = all.features,
                  OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go_DE_down, showCategory = 10, font.size = 9) + labs(title = '120 Down-regulated genes')

write.xlsx(DE, file = 'TKO_basal_vs_Ctrl_basal.xlsx', rowNames = TRUE)

#Between Ctrl and TKO differentiated populations
Idents(harmonyEpiCells) <- 'geno_cluster'
DE <- FindMarkers(harmonyEpiCells, ident.1 = c('TKO_1', 'TKO_6', 'TKO_7', 'TKO_9'), 
                  ident.2 = c('Ctrl_1', 'Ctrl_6', 'Ctrl_7'), logfc.threshold = 0.3, min.pct = 0.25)
keep <- DE$p_val_adj < 0.05 & (DE$avg_log2FC > 1 | abs(DE$pct.1 - DE$pct.2) > 0.6 | DE$avg_log2FC < -1)
topDE <- subset(DE, keep)

mat_DE <- AggrExp_Seurat@assays$RNA$scale.data[rownames(topDE),]
colnames(mat_DE) <- paste0(AggrExp_Seurat$genotype, '_', AggrExp_Seurat$Identity)
mat_DE <- mat_DE[, -c(5, 10, 19)]

Heatmap(mat_DE, name = 'TKO vs Ctrl DE genes', col = color_heatmap, 
        cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE, column_names_side = 'top', column_names_rot = 90,
        column_names_gp = gpar(col = c(rep('blue', 8), rep('red', 9))),
        row_title = '258 differentially expressed genes',
        heatmap_legend_param = list(title = 'Norm. Exp.')
)

topDE_up <- topDE[topDE$avg_log2FC > 0, ]
go_DE_up <- enrichGO(rownames(topDE_up), universe = all.features,
                     OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go_DE_up, showCategory = 10, font.size = 9) + labs(title = '187 Up-regulated genes')

topDE_down <- topDE[topDE$avg_log2FC < 0, ]
go_DE_down <- enrichGO(rownames(topDE_down), universe = all.features,
                       OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go_DE_down, showCategory = 10, font.size = 9) + labs(title = '71 Down-regulated genes')

write.xlsx(DE, file = 'TKO_diff_vs_Ctrl_diff.xlsx', rowNames = TRUE)

#Bewteen TKO hair population and Ctrl hair population
Idents(harmonyEpiCells) <- 'geno_cluster'
DE <- FindMarkers(harmonyEpiCells, ident.1 = c('TKO_3', 'TKO_5'), 
                  ident.2 = c('Ctrl_3', 'Ctrl_5'), logfc.threshold = 0.3, min.pct = 0.25)
keep <- DE$p_val_adj < 0.05 & (DE$avg_log2FC > 1 | abs(DE$pct.1 - DE$pct.2) > 0.6 | DE$avg_log2FC < -1)
topDE <- subset(DE, keep)

mat_DE <- AggrExp_Seurat@assays$RNA$scale.data[rownames(topDE),]
colnames(mat_DE) <- paste0(AggrExp_Seurat$genotype, '_', AggrExp_Seurat$Identity)
mat_DE <- mat_DE[, -c(5, 10, 19)]

Heatmap(mat_DE, name = 'TKO vs Ctrl DE genes', col = color_heatmap, 
        cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE, column_names_side = 'top', column_names_rot = 90,
        column_names_gp = gpar(col = c(rep('blue', 8), rep('red', 9))),
        row_title = '275 differentially expressed genes',
        heatmap_legend_param = list(title = 'Norm. Exp.')
)

topDE_up <- topDE[topDE$avg_log2FC > 0, ]
go_DE_up <- enrichGO(rownames(topDE_up), universe = all.features,
                     OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go_DE_up, showCategory = 10, font.size = 9) + labs(title = '127 Up-regulated genes')

topDE_down <- topDE[topDE$avg_log2FC < 0, ]
go_DE_down <- enrichGO(rownames(topDE_down), universe = all.features,
                       OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go_DE_down, showCategory = 10, font.size = 9) + labs(title = '148 Down-regulated genes')

write.xlsx(DE, file = 'TKO_diff_vs_Ctrl_diff.xlsx', rowNames = TRUE)


#A few ggplots
library(ggvenn)
RNAseq_E16_DEs <- read.xlsx('/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/RNAseq_E16_genotypeOnly_noH614_DEgenes_annotated.xlsx', 
                            sheet = 2)

DE_list <- list('BulkRNAseq' = RNAseq_E16_DEs$GeneID, 'scRNAseq' = Reduce(union, scRNA_DE))
ggvenn(DE_list, stroke_size = 1, set_name_size = 6, text_size = 5) + expand_limits(y = 2)

DE_intersect <- Reduce(intersect, DE_list)

### Between TKO_new and Ctrl Basal
Idents(harmonyEpiCells) <- 'geno_cluster'
DE <- FindMarkers(harmonyEpiCells, ident.1 = 'TKO_4', 
                  ident.2 = c('Ctrl_0', 'Ctrl_2'), logfc.threshold = 0, min.pct = 0.25)
keep <- DE$p_val_adj < 0.05 & (abs(DE$avg_log2FC) > 1 | abs(DE$pct.1 - DE$pct.2) > 0.6)
topDE <- subset(DE, keep)

mat_DE <- AggrExp_Seurat@assays$RNA$scale.data[rownames(topDE),]
colnames(mat_DE) <- paste0(AggrExp_Seurat$genotype, '_', AggrExp_Seurat$Identity)
mat_DE <- mat_DE[, -c(5, 10, 19)]

Heatmap(mat_DE, name = 'TKO vs Ctrl DE genes', col = color_heatmap, 
        cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE, column_names_side = 'top', column_names_rot = 90,
        column_names_gp = gpar(col = c(rep('blue', 8), rep('red', 9))),
        row_title = '275 differentially expressed genes',
        heatmap_legend_param = list(title = 'Norm. Exp.')
)


topDE_up <- topDE[topDE$avg_log2FC > 0, ]
go_DE_up <- enrichGO(rownames(topDE_up), universe = all.features,
                     OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go_DE_up, showCategory = 10, font.size = 9) + labs(title = '596 Up-regulated genes')

topDE_down <- topDE[topDE$avg_log2FC < 0, ]
go_DE_down <- enrichGO(rownames(topDE_down), universe = all.features,
                       OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go_DE_down, showCategory = 10, font.size = 9) + labs(title = '182 Down-regulated genes')



### Between TKO_new_2 and Ctrl Differentiated
Idents(harmonyEpiCells) <- 'geno_cluster'
DE <- FindMarkers(harmonyEpiCells, ident.1 = 'TKO_9', 
                  ident.2 = c('Ctrl_1', 'Ctrl_6', 'Ctrl_7'), logfc.threshold = 0.3, min.pct = 0.25)
keep <- DE$p_val_adj < 0.05 & (abs(DE$avg_log2FC) > 1 | abs(DE$pct.1 - DE$pct.2) > 0.6)
topDE <- subset(DE, keep)

mat_DE <- AggrExp_Seurat@assays$RNA$scale.data[rownames(topDE),]
colnames(mat_DE) <- paste0(AggrExp_Seurat$genotype, '_', AggrExp_Seurat$Identity)
mat_DE <- mat_DE[, -c(5, 10, 19)]

Heatmap(mat_DE, name = 'TKO vs Ctrl DE genes', col = color_heatmap, 
        cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE, column_names_side = 'top', column_names_rot = 90,
        column_names_gp = gpar(col = c(rep('blue', 8), rep('red', 9))),
        row_title = '275 differentially expressed genes',
        heatmap_legend_param = list(title = 'Norm. Exp.')
)

topDE_up <- topDE[topDE$avg_log2FC > 0, ]
go_DE_up <- enrichGO(rownames(topDE_up), universe = all.features,
                     OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go_DE_up, showCategory = 10, font.size = 9) + labs(title = '716 Up-regulated genes')

topDE_down <- topDE[topDE$avg_log2FC < 0, ]
go_DE_down <- enrichGO(rownames(topDE_down), universe = all.features,
                       OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go_DE_down, showCategory = 10, font.size = 9) + labs(title = '160 Down-regulated genes')


##### Keep some of the interesting DE genes for other analyses
Idents(harmonyEpiCells) <- 'geno_cluster'

DE <- FindMarkers(harmonyEpiCells, ident.1 = 'TKO_4', 
                  ident.2 = c('Ctrl_0', 'Ctrl_2'), logfc.threshold = 0.3, min.pct = 0.25)
Cluster4_DE <- subset(DE, p_val_adj < 0.05 & abs(avg_log2FC) > 1)

DE <- FindMarkers(harmonyEpiCells, ident.1 = 'TKO_9', 
                  ident.2 = c('Ctrl_1', 'Ctrl_6', 'Ctrl_7'), logfc.threshold = 0.3, min.pct = 0.25)
Cluster9_DE <- subset(DE, p_val_adj < 0.05 & abs(avg_log2FC) > 1)
Marker_C4_C9 <- union(row.names(Cluster4_DE), row.names(Cluster9_DE))

