library(DESeq2)
library(tidyverse)
library(openxlsx)

#### Performing DESeq2 ####
Counts <- Counts %>% column_to_rownames('PeakID') %>%
  select(E12_2_H3K9, E12_3_H3K9, E14_2_H3K9, E14_3_H3K9, E14_4_H3K9, E16_1_H3K9, E16_2_H3K9, E16_4_H3K9)

metaInfo <- data.frame(Name = names(Counts),
                       Stage = c('E12', 'E12', 'E14', 'E14', 'E14', 'E16', 'E16', 'E16'))
metaInfo$Stage <- factor(metaInfo$Stage, levels = c('E12', 'E14', 'E16'))
  
dds <- DESeqDataSetFromMatrix(countData = Counts, colData = metaInfo, design = ~ Stage)
peakInfo <- Epic2_union_peak
colnames(peakInfo) <- c('Chrom', 'PeakStart', 'PeakEnd', 'PeakID', 'PeakScore', 'Peak_pvalue', 'Ex_ID')
rowData(dds) <- cbind(rowData(dds), peakInfo)
rm(peakInfo)

dds <- DESeq(dds)

#### Plotting PCA
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = 'Stage')

#### Obtain dynamic peaks
Res_E16vsE12 <- results(dds, contrast = c('Stage', 'E16', 'E12'), alpha = 0.05, saveCols = c('PeakID'))
Res_E14vsE12 <- results(dds, contrast = c('Stage', 'E14', 'E12'), alpha = 0.05, saveCols = c('PeakID'))
Res_E16vsE14 <- results(dds, contrast = c('Stage', 'E16', 'E14'), alpha = 0.05, saveCols = c('PeakID'))

Dynamic_peaks <- lapply(list(Res_E16vsE12, Res_E14vsE12, Res_E16vsE14), function(res){
  df <- as.data.frame(res)
  df$padj <- ifelse(is.na(df$padj), 1, df$padj)
  dp <- subset(df, padj < 0.05 & abs(log2FoldChange) > 1)
  return(dp$PeakID)
})

Dynamic_peak_IDs <- Reduce(union, Dynamic_peaks)

#### Plot enhanced volcano plot
library(EnhancedVolcano)

df_to_plot <- as.data.frame(Res_E16vsE12)
EnhancedVolcano(df_to_plot, x = 'log2FoldChange', y = 'pvalue', lab = df_to_plot$PeakID,
                pCutoff = 0.05, pCutoffCol = 'padj', FCcutoff = 1,
                selectLab = c(''))

#### Plot the heatmap ####
library(ComplexHeatmap)

vsd <- vst(dds, blind = FALSE)
mat_z <- t(scale(t(assay(vsd))))

#color_heatmap <- colorRamp2(c(-3, 0, 3), c('blue', '#F7F7F7', 'red'))
ht <- Heatmap(mat_z[Dynamic_peak_IDs,], name = 'Dynamic Peaks',
              cluster_rows = TRUE, show_row_dend = TRUE, show_row_names = FALSE,
              cluster_columns = FALSE, column_names_side = 'top',
              row_split = 4,
              heatmap_legend_param = list(title = 'z.score'),
              use_raster = TRUE, raster_device = 'png')
ht <- draw(ht)
cluster_split <- lapply(row_order(ht), function(order){
  ids <- Dynamic_peak_IDs[order]
})

names(cluster_split) <- c('Dynamic_1', 'Dynamic_2', 'Dynamic_3', 'Dynamic_4')

#### Creating a data.frame with all DESeq comparison results and dynamic type annotations ####
FoldChangeTable_Epic2_peaks <- Epic2_union_peak %>%
  left_join(select(as.data.frame(Res_E14vsE12), log2FoldChange, padj, PeakID), by = "PeakID") %>%
  rename(FC_E14vsE12 = log2FoldChange, P_E14vsE12 = padj) %>%
  left_join(select(as.data.frame(Res_E16vsE14), log2FoldChange, padj, PeakID), by = 'PeakID') %>%
  rename(FC_E16vsE14 = log2FoldChange, P_E16vsE14 = padj) %>%
  left_join(select(as.data.frame(Res_E16vsE12), log2FoldChange, padj, PeakID), by = 'PeakID') %>%
  rename(FC_E16vsE12 = log2FoldChange, P_E16vsE12 = padj)

FoldChangeTable_Epic2_peaks <- FoldChangeTable_Epic2_peaks %>%
  mutate(Dynamics = case_when(PeakID %in% cluster_split[[1]] ~ 'Dynamic_1',
                             PeakID %in% cluster_split[[2]] ~ 'Dynamic_2',
                             PeakID %in% cluster_split[[3]] ~ 'Dynamic_3',
                             PeakID %in% cluster_split[[4]] ~ 'Dynamic_4',
                             .default = 'Constant'))

write.xlsx(FoldChangeTable_Epic2_peaks, file = 'DataFiles/FoldChangeTable_Epic2_peaks_v1.xlsx')

#### Plot heatmaps using Deeptools ####
## First obtain dynamice peaks from DESeq2 analysis
## Create bed files
write.table(subset(Epic2_union_peak, PeakID %in% cluster_split[['Dynamic_1']], select = 1:6), file = "DataFiles/Heatmap_acrossDev/DynamicPeaks_1.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
write.table(subset(Epic2_union_peak, PeakID %in% cluster_split[['Dynamic_2']], select = 1:6), file = "DataFiles/Heatmap_acrossDev/DynamicPeaks_2.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
write.table(subset(Epic2_union_peak, PeakID %in% cluster_split[['Dynamic_3']], select = 1:6), file = "DataFiles/Heatmap_acrossDev/DynamicPeaks_3.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
write.table(subset(Epic2_union_peak, PeakID %in% cluster_split[['Dynamic_4']], select = 1:6), file = "DataFiles/Heatmap_acrossDev/DynamicPeaks_4.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')


system('Deeptools --version')
cmd <- paste0('computeMatrix reference-point -S ',
              '../Bigwigs/E12_H3K9_CPM.bw ', '../Bigwigs/E14_H3K9_CPM.bw ', '../Bigwigs/E16_H3K9_CPM.bw ', '-R ',
              'DataFiles/Heatmap_acrossDev/DynamicPeaks_*.bed', ' -o ',
              'DataFiles/Heatmap_acrossDev/matrix.gz',
              ' --referencePoint center -b 5000 -a 5000 -bs 20 --sortRegions descend --smartLabels -p 20')
system(cmd)
cmd <- paste0('plotHeatmap -m ', 'DataFiles/Heatmap_acrossDev/matrix.gz', ' -o ',
              'DataFiles/Heatmap_acrossDev/Heatmap_deeptools_dynamicPeaks.pdf',
              ' --sortRegions descend --colorMap RdBu_r --whatToShow "heatmap and colorbar"',
              ' --heatmapWidth 2.5 --heatmapHeight 12'
)
system(cmd)

#### Check if up-regulated genes in bulk RNA-seq data is associated with dynamic peaks ####
RNAseq_E16_DEs <- read.xlsx('/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/RNAseq_E16_genotypeOnly_noH614_DEgenes_annotated.xlsx',
                            sheet = 1)
RNAseq_upGenes <- subset(RNAseq_E16_DEs, padj < 0.05 & log2FoldChange > 1, select = 'ensembl', drop = TRUE)

## Obtain Methylation associated genes with GO_Analysis.R script
Methyl_associated_upGenes <- lapply(Methyl_associated_genes, function(gList){
  common <- intersect(gList, RNAseq_upGenes)
})

df_to_plot <- data.frame(Dynamic = c('Constant', 'Dynamic 1', 'Dynamic 2', 'Dynamic 3', 'Dynamic 4'),
                         GeneNum = c(151, 15, 25, 6, 2)) %>%
  mutate(Percent = GeneNum / sum(GeneNum)) %>%
  mutate(Dynamic = factor(Dynamic, levels = Dynamic))

# Barplot
ggplot(df_to_plot, aes(x = Dynamic, y = GeneNum)) +
  geom_bar(stat = 'identity', fill = 'black', color = 'black', linewidth = 1, width = 0.8) +
  scale_y_log10(expand = c(0, 0.05)) +
  labs(x = NULL)

# Pie chart
ggplot(df_to_plot, aes(x = "", y = Percent, fill = Dynamic)) +
  geom_bar(stat = 'identity', width = 1) +
  coord_polar("y") +
  scale_fill_brewer(palette = 'RdYlBu', direction = -1) +
  scale_y_continuous(breaks = NULL) +
  labs(x = NULL, y = NULL, fill = 'Clusters') +
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.position = 'bottom', legend.direction = 'horizontal',
        strip.text = element_text(size = 12))


#### Check log2FoldChange of genes assocaited with different H3K9me3 dynamics
RNAseq_E16_DEs <- RNAseq_E16_DEs %>%
  mutate(Dynamics = case_when(ensembl %in% Methyl_associated_genes[['Dynamic_1']] ~ 'Dynamic_1',
                              ensembl %in% Methyl_associated_genes[['Dynamic_2']] ~ 'Dynamic_2',
                              ensembl %in% Methyl_associated_genes[['Dynamic_3']] ~ 'Dynamic_3',
                              ensembl %in% Methyl_associated_genes[['Dynamic_4']] ~ 'Dynamic_4',
                              .default = 'Constant'))
ggplot(RNAseq_E16_DEs, aes(x = Dynamics, y = log2FoldChange)) +
  geom_boxplot()


