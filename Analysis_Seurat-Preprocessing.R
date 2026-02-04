library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
#library(xlsx)

folder_path <- "/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/scRNA-seq/CellRangerOutPuts/"

Data_CT1 <- Read10X(paste0(folder_path, "CT_1/raw_feature_bc_matrix"))
Data_CT2 <- Read10X(paste0(folder_path, "CT_2/raw_feature_bc_matrix"))
Data_KO1 <- Read10X(paste0(folder_path, "KO_1/raw_feature_bc_matrix"))
Data_KO2 <- Read10X(paste0(folder_path, "KO_2/raw_feature_bc_matrix"))

raw_CT_1 <- CreateSeuratObject(counts = Data_CT1, min.cells = 10, min.features = 500)
raw_CT_2 <- CreateSeuratObject(counts = Data_CT2, min.cells = 10, min.features = 500)
raw_KO_1 <- CreateSeuratObject(counts = Data_KO1, min.cells = 10, min.features = 500)
raw_KO_2 <- CreateSeuratObject(counts = Data_KO2, min.cells = 10, min.features = 500)

rm(list = ls()[grep("^Data", ls())])

allCells <- merge(raw_CT_1, y = c(raw_CT_2, raw_KO_1, raw_KO_2), add.cell.id = c("CT1", "CT2", "KO1", "KO2"), project = "E16_basalCells")
allCells$orig.ident <- c(rep("CT1", 5517), rep("CT2", 5218), rep("KO1", 2951), rep("KO2", 3940))
allCells$genotype <- c(rep("Ctrl", 5517), rep("Ctrl", 5218), rep("TKO", 2951), rep("TKO", 3940))

rm(list = ls()[grep("^raw", ls())])

#First glance and filter
allCells[["percent.mt"]] <- PercentageFeatureSet(allCells, pattern = "^mt-")
VlnPlot(allCells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), split.by = "genotype", group.by = "orig.ident",
        alpha = 0, cols = c('darkgreen', 'gold'))
FeatureScatter(allCells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#Data pre-processing
filteredCells <- subset(allCells, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)
filteredCells <- NormalizeData(filteredCells, normalization.method = "LogNormalize", scale.factor = 10000)
filteredCells <- FindVariableFeatures(filteredCells, selection.method = "vst", nfeatures = 2000)
all.features <- rownames(filteredCells)
filteredCells <- ScaleData(filteredCells, features = all.features)
VlnPlot(filteredCells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), split.by = "genotype", group.by = "orig.ident",
        alpha = 0, cols = c('darkgreen', 'gold'))

#Dimentional reduction and cluster analysis (without integration)
filteredCells <- RunPCA(filteredCells, features = VariableFeatures(filteredCells), reduction.name = 'pca')
ElbowPlot(filteredCells)
filteredCells <- FindNeighbors(filteredCells, dims = 1:20)

filteredCells <- FindClusters(filteredCells, resolution = 0.25)
filteredCells <- RunUMAP(filteredCells, dims = 1:20, reduction.name = 'umap_no_int')
DimPlot(filteredCells, reduction = "umap_no_int", split.by = "genotype", pt.size = 0.6, label = TRUE)

cluster_count <- data.frame(sample = filteredCells$genotype, cluster = Idents(filteredCells)) %>% 
  group_by(sample, cluster) %>% summarise(count = n(), .groups = 'drop')

#Cell markers
marker_immune <- c('Spi1', 'Cd86', 'Cx3cr1', 'Ptprc')
marker_fibroblast <- c('Col1a1', 'Pdgfra', 'Fbn1', 'Twist2')
marker_melanocyte <- c('Dct', 'Tyr', 'Kit', 'Mitf')
marker_merkel <- c('Krt8', 'Sox2')
marker_keratinocyte <- c('Krt14', 'Krt10', 'Trp63', 'Epcam')
marker_cluster9 <- c('Sox2', 'Prox1', 'Syt7', 'Krt8', 'Map2', 'Piezo2')

FeaturePlot(filteredCells, features = c('Krt14', 'Krt10'), pt.size = 0.7, cols = c('lightgrey', 'red'))

#Initial assignment of cluster identity
DotPlot(filteredCells, features = c(marker_immune, marker_fibroblast, marker_melanocyte, marker_cluster9), 
        dot.scale = 6, col.min = -5, col.max = 5) +
  RotatedAxis()

filteredCells_layersJoint <- filteredCells
filteredCells_layersJoint[['RNA']] <- JoinLayers(filteredCells[['RNA']])

marker_9 <- FindConservedMarkers(filteredCells_layersJoint, ident.1 = '9', grouping.var = 'genotype', min.cells.group = 10)
marker_10 <- FindConservedMarkers(filteredCells_layersJoint, ident.1 = '10', grouping.var = 'genotype', min.cells.group = 10)
marker_9 <- subset(marker_9, marker_9$Ctrl_p_val_adj < 0.05 & marker_9$Ctrl_avg_log2FC > 2)
marker_10 <- subset(marker_10, marker_10$TKO_p_val_adj < 0.05 & marker_10$TKO_avg_log2FC > 2)

rm(filteredCells_layersJoint)

#Cluster 5 and 11 are immune population;
DimPlot(filteredCells, reduction = 'umap_no_int', pt.size = 0.7,
        cells.highlight = list('5'=WhichCells(filteredCells, idents = 5), '11'=WhichCells(filteredCells, idents = 11)),
        cols.highlight = c('darkgreen', 'lightgreen'), split.by = 'genotype', ncol = 1)
FeaturePlot(filteredCells, features = marker_immune, pt.size = 0.7, cols = c('lightgrey', 'red'))
VlnPlot(filteredCells, features = marker_immune, idents = c('5', '11'), split.by = 'genotype', ncol = 2)
#Cluster 6 is fibroblast;
DimPlot(filteredCells, reduction = 'umap_no_int', pt.size = 0.7,
        cells.highlight = list('6'=WhichCells(filteredCells, idents = 6)),
        cols.highlight = 'blue', split.by = 'genotype', ncol = 1)
FeaturePlot(filteredCells, features = marker_fibroblast, pt.size = 0.7, cols = c('lightgrey', 'red'))
#Cluster 7 is TKO specific
#Cluster 8 is diminished in TKO
DimPlot(filteredCells, reduction = 'umap_no_int', pt.size = 0.7,
        cells.highlight = list('7'=WhichCells(filteredCells, idents = 7), '8'=WhichCells(filteredCells, idents = 8)),
        cols.highlight = c('darkgreen', 'blue'), split.by = 'genotype', ncol = 1)
VlnPlot(filteredCells, features = c('Krt14', 'Trp63', 'Krt10', 'Itga6'), idents = c('7', '8'), split.by = 'genotype', ncol = 2)

#Cluster 9 is control specific, presumably merkel cell population
DimPlot(filteredCells, reduction = 'umap_no_int', pt.size = 0.7,
        cells.highlight = list('9'=WhichCells(filteredCells, idents = 9)),
        cols.highlight = 'brown', split.by = 'genotype', ncol = 1)
FeaturePlot(filteredCells, features = marker_cluster9, pt.size = 0.7, cols = c('lightgrey', 'red'), ncol = 3)
FeaturePlot(filteredCells, features = marker_merkel, pt.size = 0.7, cols = c('lightgrey', 'red'), ncol = 3)
go_marker_9 <- enrichGO(rownames(marker_9), OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go_marker_9, showCategory = 10)

#Cluster 10 is melanocyte population
DimPlot(filteredCells, reduction = 'umap_no_int', pt.size = 0.7,
        cells.highlight = list('10'=WhichCells(filteredCells, idents = 10)),
        cols.highlight = 'violet', split.by = 'genotype', ncol = 1)
FeaturePlot(filteredCells, features = marker_melanocyte, pt.size = 0.7, cols = c('lightgrey', 'red'))
go_marker_10 <- enrichGO(rownames(marker_10), OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go_marker_10, showCategory = 10)

#Assigning identities
Ident_no_int_1st <- c('0', '1', '2', '3', '4', 'Cd86+ Immune', 'Fibroblasts', '7', '8', 'Merkel cells', 'Melanocytes', 'Cd86- Immune')
names(Ident_no_int_1st) <- levels(filteredCells)
filteredCells <- RenameIdents(filteredCells, Ident_no_int_1st)
filteredCells$Ident_no_int <- Idents(filteredCells)

DimPlot(filteredCells, reduction = "umap_no_int", pt.size = 3, split.by = 'genotype', label = TRUE,
        raster = TRUE)
DotPlot(filteredCells, dot.min = 0.2,
        features = c(marker_keratinocyte, marker_immune, marker_fibroblast, marker_melanocyte, marker_merkel),
        cluster.idents = TRUE) +
  RotatedAxis()

#Removing unwanted populations
filteredCells_named$clusters_no_int_res0.25 <- Idents(filteredCells_named)
EpidermalCells <- subset(filteredCells_named, idents = c('0', '1', '2', '3', '4', '7', '8'))

saveRDS(filteredCells, file = "Filtered cells.RData")
saveRDS(EpidermalCells, file = "Epidermal cells.RData")

rm(filteredCells_named)
rm(filteredCells)
rm(EpidermalCells)

save.image(file = "Preprocessing.RData")

#Loading saved seurat object for re-work
filteredCells <- readRDS(file = "Filtered cells.RData")
EpidermalCells <- readRDS(file = "Epidermal cells.RData")


############################################################################
#Import data generated by cellranger_aggr
scRNA_Data <- Read10X("C:/Users/yhbai/OneDrive - Cornell University/Tumbar Lab/scRNA-seq/CellRangerOutPuts/count/filtered_feature_bc_matrix")
allCells <- CreateSeuratObject(counts = scRNA_Data, min.cells = 10, min.features = 500)