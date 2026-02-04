library(Seurat)
library(tidyverse)
library(patchwork)
library(openxlsx)

saveRDS(harmonyEpiCells, file = "Harmony_integrated_epidermal_w_Merkel.RData")
harmonyEpiCells <- readRDS(file = "Harmony_integrated_epidermal_w_Merkel.RData")


## Remove small clusters (cell number < 10), which is unlikely real and creates problems in downstream analyses
Idents(harmonyEpiCells) <- 'geno_cluster'
harmonyEpiCells <- subset(harmonyEpiCells, ident = c('Ctrl_4', 'Ctrl_9', 'TKO_8'), invert = TRUE)

## Cluster definition
Idents(harmonyEpiCells) <- 'seurat_clusters'
Cluster_idents <- c('0' = 'Basal_1', '1' = 'Diff_1', '2' = 'Basal_2', '3' = 'HairFollicle',
                    '4' = 'AberrBasal', '5' = 'Basal_3', '6' = 'Diff_2', '7' = 'Diff_3',
                    '8' = 'Merkel', '9' = 'AberrDiff')
harmonyEpiCells <- RenameIdents(harmonyEpiCells, Cluster_idents)

harmonyEpiCells$genotype <- factor(harmonyEpiCells$genotype, levels = c('Ctrl', 'TKO'))
harmonyEpiCells$Identity <- factor(Idents(harmonyEpiCells),
                                   levels = c('Basal_1', 'Basal_2', 'Basal_3', "AberrBasal",
                                              'Diff_1', 'Diff_2', 'Diff_3', 'AberrDiff',
                                              'HairFollicle', 'Merkel'))
Idents(harmonyEpiCells) <- 'Identity'

cluster_colors <- rev(RColorBrewer::brewer.pal(10, 'Set3'))
cluster_colors[2] <- rgb(217, 187, 140, maxColorValue = 255)
cluster_colors[3] <- rgb(252, 214, 17, maxColorValue = 255)

DimPlot(harmonyEpiCells, reduction = 'umap', split.by = 'genotype', label.size = 5,
        label = FALSE, raster = FALSE) +
  scale_color_manual(values = rev(cluster_colors))

all.features <- rownames(harmonyEpiCells)

#Dotplot showing marker genes validating cell identities
marker_genes <- c("Trp63", "Krt14", "Itga6", "Mki67", "Cdk1", "Krt10", "Krt1", "Elovl4",
                  "Gli2", "Eda", "Lhx2", "Sox9", "Runx1", "Krt8", "Krt20")

plot <- DotPlot(harmonyEpiCells, features = marker_genes, cols = 'RdBu') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7),
        axis.title = element_blank())

plot$data$id <- factor(plot$data$id, levels = rev(cluster_levels))

#### PART I: characterize population-level changes (Appearance of two aberrant population) ####
wb <- loadWorkbook("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/scRNA_DE_genes_v4_annotated.xlsx")
ClusterDEgenes <- list()
for (name in names(wb)) {
  ClusterDEgenes[[name]] <- readWorkbook(wb, sheet = name)
}
rm(wb, name)

library(clusterProfiler)
library(org.Mm.eg.db)

geneList <- ClusterDEgenes[['4']]
genes_up <- unique(subset(geneList, p_val_adj < 0.05 & avg_log2FC > 1, select = GeneID, drop = TRUE))
genes_down <- unique(subset(geneList, p_val_adj < 0.05 & avg_log2FC < -1, select = GeneID, drop = TRUE))

go_up <- enrichGO(genes_up, universe = all.features, OrgDb = org.Mm.eg.db,
                  ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
go_down <- enrichGO(genes_down, universe = all.features, OrgDb = org.Mm.eg.db,
                    ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
# Write GO result to excel files for manual picking
local({
  wb <- createWorkbook()
  addWorksheet(wb, "Genes_Up")
  addWorksheet(wb, "Genes_Down")
  writeData(wb, sheet = "Genes_Up", as.data.frame(go_up), rowNames = TRUE)
  writeData(wb, sheet = "Genes_Down", as.data.frame(go_down), rowNames = TRUE)
  saveWorkbook(wb, file = "ExcelFiles/GO_Markers_AberrBasal.xlsx", overwrite = TRUE)
})

geneList <- ClusterDEgenes[['9']]
genes_up <- unique(subset(geneList, p_val_adj < 0.05 & avg_log2FC > 1, select = GeneID, drop = TRUE))
genes_down <- unique(subset(geneList, p_val_adj < 0.05 & avg_log2FC < -1, select = GeneID, drop = TRUE))

go_up <- enrichGO(genes_up, universe = all.features, OrgDb = org.Mm.eg.db,
                  ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
go_down <- enrichGO(genes_down, universe = all.features, OrgDb = org.Mm.eg.db,
                    ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
# Write GO result to excel files for manual picking
local({
  wb <- createWorkbook()
  addWorksheet(wb, "Genes_Up")
  addWorksheet(wb, "Genes_Down")
  writeData(wb, sheet = "Genes_Up", as.data.frame(go_up), rowNames = TRUE)
  writeData(wb, sheet = "Genes_Down", as.data.frame(go_down), rowNames = TRUE)
  saveWorkbook(wb, file = "ExcelFiles/GO_Markers_AberrDiff.xlsx", overwrite = TRUE)
})

rm(go_up, go_down, geneList, genes_up, genes_down)

## Plot manually picked go terms ##
GO_terms_to_plot <- c("GO:0007264", "GO:0031032", "GO:0031589", "GO:0006007", "GO:0007015",
                      "GO:0016358", "GO:0034329", "GO:0090559", "GO:0099173", "GO:0003013")
go_table <- read.xlsx("ExcelFiles/GO_Markers_AberrBasal.xlsx", sheet = 1, rowNames = TRUE)
p1 <- CB_plot_GO(GO_terms = GO_terms_to_plot, go_table = go_table) +
  labs(title = "Upregulated") +
  theme(axis.title.x = element_blank(),
        legend.position = 'none')

GO_terms_to_plot <- c("GO:0008544", "GO:0043588", "GO:0030198", "GO:0009611", "GO:0006909",
                      "GO:0002685", "GO:0061564", "GO:0021782", "GO:1905517", "GO:0014002")
go_table <- read.xlsx("ExcelFiles/GO_Markers_AberrBasal.xlsx", sheet = 2, rowNames = TRUE)
p2 <- CB_plot_GO(GO_terms = GO_terms_to_plot, go_table = go_table) +
  labs(title = "Downregulated")

GO_terms_to_plot <- c("GO:0031589", "GO:0010631", "GO:0001558", "GO:0009611", "GO:0097191",
                      "GO:0016049", "GO:0003013", "GO:0097421", "GO:0043588", "GO:0016055")
go_table <- read.xlsx("ExcelFiles/GO_Markers_AberrDiff.xlsx", sheet = 1, rowNames = TRUE)
p3 <- CB_plot_GO(GO_terms = GO_terms_to_plot, go_table = go_table) +
  labs(title = "Upregulated") +
  theme(axis.title.x = element_blank(),
        legend.position = 'none')

GO_terms_to_plot <- c("GO:0030216", "GO:0009913", "GO:0006631", "GO:0043588", "GO:0006633",
                      "GO:0046394", "GO:0052547", "GO:0046512", "GO:0050727", "GO:0006641")
go_table <- read.xlsx("ExcelFiles/GO_Markers_AberrDiff.xlsx", sheet = 2, rowNames = TRUE)
p4 <- CB_plot_GO(GO_terms = GO_terms_to_plot, go_table = go_table) +
  labs(title = "Downregulated")

p1 + p2 + p3 + p4

rm(p1, p2, p3, p4, GO_terms_to_plot, go_table)

## Use simplifyEnrichment package
library(simplifyEnrichment)
#Upregulated genes
go_table <- list('AberrBasal_Up' = read.xlsx('ExcelFiles/GO_Markers_AberrBasal.xlsx', sheet = 1, rowNames = TRUE),
                 'AberrDiff_Up' = read.xlsx('ExcelFiles/GO_Markers_AberrDiff.xlsx', sheet = 1, rowNames = TRUE),
                 'AberrBasal_Down' = read.xlsx('ExcelFiles/GO_Markers_AberrBasal.xlsx', sheet = 2, rowNames = TRUE),
                 'AberrDiff_Down' = read.xlsx('ExcelFiles/GO_Markers_AberrDiff.xlsx', sheet = 2, rowNames = TRUE))

go_table <- lapply(go_table, function(dt){subset(dt, ONTOLOGY =='BP')})

pdf("Plots/GO_Clusters_AberrPop.pdf", width = 9, height = 6)
simplifyGOFromMultipleLists(lt = go_table, padj_cutoff = 0.05)
dev.off()


#### PART II: show TKO-induced transcriptomic changes are cluster specific

# Prepare expression matrix
AggrExp_Seurat <- AggregateExpression(harmonyEpiCells, return.seurat = TRUE, group.by = c('genotype', 'seurat_clusters'))
AggrExp_Seurat$geno_cluster <- Idents(AggrExp_Seurat)
Idents(AggrExp_Seurat) <- 'seurat_clusters'
AggrExp_Seurat <- RenameIdents(AggrExp_Seurat, Cluster_idents)
AggrExp_Seurat$Identity <- Idents(AggrExp_Seurat)

mat_DE <- AggrExp_Seurat@assays$RNA$scale.data

#Plot RNA-seq DE genes
RNAseq_E16_DEs <- read.xlsx('/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/RNAseq_E16_genotypeOnly_noH614_DEgenes_annotated.xlsx', 
                            sheet = 2)

features_to_plot <- RNAseq_E16_DEs %>% select(c(GeneID, log2FoldChange)) %>% arrange(desc(log2FoldChange)) %>%
  mutate(Reg = if_else(log2FoldChange > 1, 'Upregulated', 'Downregulated')) %>% subset(GeneID %in% all.features)

## Make the heatmap
library(ComplexHeatmap)
library(circlize)

color_heatmap <- colorRamp2(c(-3, 0, 3), c('blue', '#F7F7F7', 'red'))

Heatmap(mat_DE[features_to_plot$GeneID, ], name = 'TKO vs Ctrl DE genes', col = color_heatmap, 
        cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE, row_names_side = 'left',
        column_order = c(13,17,15,16,10,14,12,11,9,1,3,5,4,2,7,6,8),
        row_split = features_to_plot$Reg, row_title_rot = 0,
        column_names_side = 'top', column_names_rot = 90, column_labels = AggrExp_Seurat$Identity,
        column_names_gp = gpar(col = c(rep('blue', 8), rep('red', 9))),
        #row_title = 'Differentially expressed genes',
        heatmap_legend_param = list(title = 'Norm. Exp.'),
        use_raster = TRUE, raster_device = 'png'
)

# Plot top markers for each cluster
wb <- loadWorkbook("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/scRNA_Conserved_markers_v3_annotated.xlsx")
ClusterMarkers <- list()
for (name in names(wb)) {
  ClusterMarkers[[name]] <- readWorkbook(wb, sheet = name, rowNames = TRUE)
}
rm(wb, name)

topMarkers <- lapply(ClusterMarkers, function(dt){
  names <- rownames(dt)
  if ("TKO_avg_log2FC" %in% colnames(dt)){
    top <- dt %>% subset(TKO_p_val_adj < 0.05) %>% top_n(50, TKO_avg_log2FC)
    names <- intersect(names, rownames(top))
  }
  if ("Ctrl_avg_log2FC" %in% colnames(dt)){
    top <- dt %>% subset(Ctrl_p_val_adj < 0.05) %>% top_n(50, Ctrl_avg_log2FC)
    names <- intersect(names, rownames(top))
  }
  return(names)
})
features_to_plot <- unique(Reduce(union, topMarkers))

genes_of_interest <- c('Krt14', 'Itga6', 'Mki67', 'Cdk1', 'Krt10', 'Flnb', 'Sox9', 'Lhx2', 'Krt8')
features_to_plot <- unique(union(features_to_plot, genes_of_interest))
ha <- rowAnnotation(geneLabel = anno_mark(at = match(genes_of_interest, features_to_plot), labels = genes_of_interest))
ht <-
  Heatmap(mat_DE[features_to_plot, ], name = 'TKO vs Ctrl DE genes', col = color_heatmap, right_annotation = ha,
          cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE, row_names_side = 'left',
          column_order = c(1,9,3,11,5,14,4,12,2,10,6,15,7,16,13,8,17),
          column_names_side = 'top', column_names_rot = 90, column_labels = AggrExp_Seurat$Identity,
          column_names_gp = gpar(col = c(rep('blue', 8), rep('red', 9))),
          row_title = 'Marker Genes',
          heatmap_legend_param = list(title = 'Norm. Exp.', direction = 'horizontal'),
          use_raster = TRUE, raster_device = 'png'
 )
draw(ht, heatmap_legend_side = 'bottom')
rm(ha, ht)

# Plot transposed heatmap
Heatmap(t(mat_DE[features_to_plot, ]), name = 'TKO vs Ctrl DE genes', col = color_heatmap, 
        cluster_columns = FALSE, show_column_dend = FALSE, show_column_names = FALSE, column_names_side = 'top', column_km = 11, cluster_column_slices = FALSE,
        row_order = c(1,9,3,11,5,14,4,12,2,10,6,15,7,16,13,8,17),
        row_names_side = 'left', row_labels = AggrExp_Seurat$Identity,
        row_names_gp = gpar(col = c(rep('blue', 8), rep('red', 9))),
        column_title = 'Marker Genes',
        heatmap_legend_param = list(title = 'Norm. Exp.'),
        use_raster = TRUE, raster_device = 'png'
)


# Plot top DE genes obtained from single cell anlaysis
scRNA_DE <- lapply(ClusterDEgenes[-c(5,9,10)], function(dt){
  names <- subset(dt, p_val_adj < 0.05 & abs(avg_log2FC) > 1) %>% select(c(GeneID, avg_log2FC)) %>%
    arrange(desc(avg_log2FC)) %>% mutate(Reg = if_else(avg_log2FC > 0, 'Upregulated', 'Downregulated'))
})

features_to_plot <- bind_rows(scRNA_DE) %>% distinct(GeneID, .keep_all = TRUE)
Heatmap(mat_DE[features_to_plot$GeneID, ], name = 'TKO vs Ctrl DE genes', col = color_heatmap, 
        cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE, row_names_side = 'left',
        #column_order = c(1,9,3,11,5,14,4,12,2,10,6,15,7,16,13,8,17),
        row_split = features_to_plot$Reg, row_title_rot = 90, row_title_gp = gpar(fontsize = 12),
        column_names_side = 'top', column_names_rot = 90, column_labels = AggrExp_Seurat$Identity,
        column_names_gp = gpar(col = c(rep('blue', 8), rep('red', 9))),
        heatmap_legend_param = list(title = 'Norm. Exp.'),
        use_raster = TRUE, raster_device = 'png'
)

# Plot DE obtained by comparing TKO vs Ctrl in Seurat
Idents(harmonyEpiCells) <- 'genotype'
DE <- FindMarkers(harmonyEpiCells, ident.1 = 'TKO', ident.2 = 'Ctrl', logfc.threshold = 0.3, min.pct = 0.25)
features_to_plot <- rownames(subset(DE, abs(avg_log2FC) > 1 & p_val_adj < 0.05))

Heatmap(mat_DE[features_to_plot, ], name = 'TKO vs Ctrl DE genes', col = color_heatmap, 
        cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE, row_names_side = 'left',
        #column_order = c(1,9,3,11,5,14,4,12,2,10,6,15,7,16,13,8,17),
        column_names_side = 'top', column_names_rot = 90, column_labels = AggrExp_Seurat$Identity,
        column_names_gp = gpar(col = c(rep('blue', 8), rep('red', 9))),
        row_title = 'Differentially expressed genes',
        heatmap_legend_param = list(title = 'Norm. Exp.'),
        use_raster = TRUE, raster_device = 'png'
)


### Plot DE genes from from Nascent Analysis
Nascent_fullGene_DE <- read.xlsx('/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/Precise2024_DiffExpressedGenes_fullGene.xlsx', rowNames = TRUE)
Nascent_up_genes <- subset(Nascent_fullGene_DE, padj < 0.05 & log2FoldChange > 1, select = 'SYMBOL', drop = TRUE)

features_to_plot <- intersect(Nascent_up_genes, all.features)
Heatmap(mat_DE[features_to_plot, ], name = 'TKO vs Ctrl DE genes', col = color_heatmap, 
        cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE, row_names_side = 'left',
        #column_order = c(1,9,3,11,5,14,4,12,2,10,6,15,7,16,13,8,17),
        column_names_side = 'top', column_names_rot = 90, column_labels = AggrExp_Seurat$Identity,
        column_names_gp = gpar(col = c(rep('blue', 8), rep('red', 9))),
        row_title = 'Differentially expressed genes',
        heatmap_legend_param = list(title = 'Norm. Exp.'),
        use_raster = TRUE, raster_device = 'png'
)


### Plot marker genes that show mixed lineage identity
genes_aberrDiff <- as.data.frame(mat_DE) %>% top_n(5000, TKO_9)

mat_to_plot <-mat_DE[rownames(genes_aberrDiff), c(1:8, 17)]

Heatmap(mat_to_plot, name = 'TKO vs Ctrl DE genes', col = color_heatmap, 
        cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE, row_names_side = 'left',
        #column_order = c(1,9,3,11,5,14,4,12,2,10,6,15,7,16,13,8,17),
        cluster_columns = TRUE, column_names_side = 'top', column_names_rot = 90, column_labels = AggrExp_Seurat$Identity[c(1:8,17)],
        column_names_gp = gpar(col = c(rep('blue', 8), rep('red', 9))),
        row_title = 'Differentially expressed genes',
        heatmap_legend_param = list(title = 'Norm. Exp.'),
        use_raster = TRUE, raster_device = 'png'
)

rm(mat_to_plot)

## Subset Seurat object to find marker genes to demonstrate mixed lineage identity
Idents(harmonyEpiCells) <- 'geno_cluster'
Sub_Seurat <- subset(harmonyEpiCells, ident = c('Ctrl_0', 'Ctrl_1', 'Ctrl_2', 'Ctrl_3', 'Ctrl_5',
                                                'Ctrl_6', 'Ctrl_7', 'Ctrl_8', 'TKO_9'), invert = FALSE)
Idents(harmonyEpiCells) <- 'Identity'

Marker_list <- list()
for (cluster in c('Ctrl_0', 'Ctrl_1', 'Ctrl_2', 'Ctrl_3', 'Ctrl_5', 'Ctrl_6', 'Ctrl_7', 'Ctrl_8')){
  marker <- FindMarkers(Sub_Seurat, ident.1 = c(cluster, 'TKO_9'), logfc.threshold = 0.3, min.pct = 0.25)
  Marker_list[[cluster]] <- marker
}
rm(cluster, marker)

topMarker_list <- lapply(Marker_list, function(marker){
  df <- marker %>% subset(p_val_adj < 0.05) %>% top_n(5000, avg_log2FC)
  name <- intersect(rownames(df), rownames(genes_aberrDiff))
})
topMarker_list[['TKO_9']] <- ClusterMarkers[['9']] %>% subset(TKO_p_val_adj < 0.05) %>% top_n(500, TKO_avg_log2FC) %>% rownames()

features_to_plot <- unique(Reduce(union, topMarker_list))
mat_to_plot <- mat_DE[features_to_plot, c(1:8, 13, 17)]

Heatmap(mat_to_plot, name = 'TKO vs Ctrl DE genes', col = color_heatmap, 
        cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE, row_names_side = 'left',
        #column_order = c(1,9,3,11,5,14,4,12,2,10,6,15,7,16,13,8,17),
        cluster_columns = TRUE, column_names_side = 'top', column_names_rot = 90, column_labels = AggrExp_Seurat$Identity[c(1:8, 13, 17)],
        column_split = c(rep('Ctrl', 8), rep('TKO', 2)),
        column_names_gp = gpar(col = c(rep('blue', 8), rep('red', 9))),
        row_title = 'AberrDiff marker genes',
        heatmap_legend_param = list(title = 'Norm. Exp.'),
        use_raster = TRUE, raster_device = 'png'
)

#Plot horizontally
Heatmap(t(mat_to_plot), name = 'TKO vs Ctrl DE genes', col = color_heatmap, 
        cluster_columns = TRUE, show_column_dend = FALSE, show_column_names = FALSE, column_names_side = 'bottom',
        #column_order = c(1,9,3,11,5,14,4,12,2,10,6,15,7,16,13,8,17),
        cluster_rows = TRUE, row_names_side = 'left', column_names_rot = 0, row_labels = AggrExp_Seurat$Identity[c(1:8, 13, 17)],
        row_split = factor(c(rep('Ctrl', 8), rep('TKO', 2)), levels = c('TKO', 'Ctrl')),
        row_names_gp = gpar(col = c(rep('blue', 8), rep('red', 2))),
        column_title = 'AberrDiff marker genes',
        heatmap_legend_param = list(title = 'Norm. Exp.'),
        use_raster = TRUE, raster_device = 'png'
)

rm(topMarker_list, Marker_list)
rm(Sub_Seurat)


