###This is primarily additional heatmap plotting, if missing anything, trying finding
library(tidyverse)
library(openxlsx)

#Loading some of the necessary files
txdb_ensembl <- makeTxDbFromEnsembl(organism = 'Mus musculus', release = 110)
all_genes <- genes(txdb_ensembl)
info <- select(Ens_mm39_110, keys = unlist(all_genes$gene_id), keytype = 'GENEID', 
               columns = c('GENEBIOTYPE', 'GENEID', 'GENENAME', 'SYMBOL', 'DESCRIPTION'))
all_genes$Biotype <- info$GENEBIOTYPE
all_genes$Symbol <- info$SYMBOL
rm(info)
seqlevelsStyle(all_genes) <- 'UCSC'
all_genes <- keepSeqlevels(all_genes, chrs, pruning.mode = 'tidy')
all_genes <- sort(all_genes)
names(all_genes) <- NULL

Methyl_Clusters <- read.xlsx("DataFiles/Methyl_Annotation/Methyl_clusters.xlsx", colNames = TRUE)

system2("deeptools", args = "--version")

####### Ploting for upregulated genes from different single cell clusters ######
wb <- loadWorkbook("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/scRNA_Conserved_markers_v3_annotated.xlsx")
ClusterMarkers <- list()
for (name in names(wb)) {
  ClusterMarkers[[name]] <- readWorkbook(wb, sheet = name)
}

wb <- loadWorkbook("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/scRNA_DE_genes_v4_annotated.xlsx")
ClusterDEgenes <- list()
for (name in names(wb)) {
  ClusterDEgenes[[name]] <- readWorkbook(wb, sheet = name)
}

rm(wb, name)

Cluster_idents <- c('0' = 'Basal_1', '1' = 'Diff_1', '2' = 'Basal_2',
                    '3' = 'Hair_follicle', '4' = 'AberrBasal', '5' = 'Basal_3',
                    '6' = 'Diff_2', '7' = 'Diff_3', '8' = 'Merkel',
                    '9' = 'AberrDiff')

DEgenes_Up <- lapply(ClusterDEgenes, function(df) df[df$avg_log2FC > 0,1])
DEgenes_Down <- lapply(ClusterDEgenes, function(df) df[df$avg_log2FC < 0,1])
Up_genes_combined <- Reduce(union, DEgenes_Up[-c(5,9,10)])
Down_genes_combined <- Reduce(union, DEgenes_Down[-c(5,9,10)])

Markers_Up <- lapply(ClusterMarkers[c('4','9')], function(df) df[df$TKO_avg_log2FC > 3,1])
Markers_Down <- lapply(ClusterMarkers[c('4','9')], function(df) df[df$TKO_avg_log2FC < -3,1])
TKO_Markers_Up <- Reduce(union, Markers_Up)
TKO_Markers_Down <- Reduce(union, Markers_Down)

scRNA_GOI <- DEgenes_Up
names(scRNA_GOI) <- Cluster_idents[names(scRNA_GOI)]

Cluster_GOI <- list()
for (name in names(scRNA_GOI)) {
  gList <- scRNA_GOI[[name]]
  Cluster_GOI[[name]] <- Methyl_Clusters %>% subset(subset = Ensembl %in% gList | Symbol %in% gList, select = c(Cluster, Symbol, Ensembl)) %>%
    mutate(CellCluster = name)
}
rm(name, gList)

df_to_plot <- do.call(rbind, Cluster_GOI) %>% group_by(CellCluster, Cluster) %>% summarise(count = n()) %>%
  mutate(Percent = count / sum(count)) %>% ungroup
df_to_plot <- Methyl_Clusters %>% group_by(Cluster) %>% summarise(count = n()) %>% 
  mutate(Percent = count / sum(count), CellCluster = 'AllGenes') %>% ungroup() %>%
  rbind(df_to_plot, .) %>%
  mutate(Percent = if_else(Cluster == 'C5', NA, Percent))

ggplot(df_to_plot, aes(x = CellCluster, y = Percent, fill = Cluster)) +
  geom_bar(stat = 'identity', position = position_stack(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0.005)) +
  scale_fill_brewer(palette = 'RdYlBu')+
  labs(y = '% Genes', fill = 'Cluster', x = element_blank()) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.8))

for (name in names(Cluster_GOI)) {
  cat("Processing ", name, "\n")
  CB_pHeatmap_methylCluster(name, path = "DataFiles/ClusterHeatmaps/scRNAupGenes/")
}


####### Ploting for different gene categories from nascent transcription data #####
PreciseSeqFC_table <- read.xlsx("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/Precise2024_FoldChangeTable_v2.xlsx")
PreciseSeqFC_table <- PreciseSeqFC_table %>% mutate(Category = factor(Category,
                                      levels = c('upPI_nsGB', 'nsPI_upGB', 'downPI_upGB', 'downPI_nsGB', 'others')))
Nascent_GOI <- list()
for (name in levels(PreciseSeqFC_table$Category)){
  Nascent_GOI[[name]] <- subset(PreciseSeqFC_table, Category == name, select = 'GeneID', drop = TRUE)
}

for (name in names(Nascent_GOI)){
  gr <- subset(all_genes, gene_id %in% Nascent_GOI[[name]])
  CB_export_bed(gr, paste0("DataFiles/Methyl_Plotting/Genebody_Set_", name, ".bed"), identifier = 'gene_id')
}

#Plot without clustering
command <- paste0("computeMatrix reference-point -S ", "DataFiles/H3K9me3_IgG_subtracted.bw",
                  " -R ",
                  "DataFiles/Methyl_Plotting/Genebody_Set_upPI_nsGB.bed ", "DataFiles/Methyl_Plotting/Genebody_Set_nsPI_upGB.bed ",
                  "DataFiles/Methyl_Plotting/Genebody_Set_downPI_upGB.bed ", "DataFiles/Methyl_Plotting/Genebody_Set_downPI_nsGB.bed ",
                  "DataFiles/Methyl_Plotting/Genebody_Random.bed ",
                  " -o ", "DataFiles/Methyl_Plotting/matrix_nascent_set.gz",
                  " --referencePoint TSS --binSize 200 --sortRegions descend -p 15 -b 10000 -a 10000")
system(command)
command <- paste("plotHeatmap -m ", "DataFiles/Methyl_Plotting/matrix_nascent_set.gz",
                 " -o ", "DataFiles/Methyl_Plotting/Methylation_aroundTSS_Nascent_set.pdf",
                 " --colorMap RdBu_r --samplesLabel H3K9me3 --regionsLabel upPI_nsGB nsPI_upGB downPI_upGB downPI_nsGB random",
                 " --zMin -0.25 --zMax 0.25 --heatmapWidth 8 --heatmapHeight 16")
system(command)

#Plot with clustering
Cluster_GOI <- list()
for (name in names(Nascent_GOI)) {
  gList <- Nascent_GOI[[name]]
  Cluster_GOI[[name]] <- Methyl_Clusters %>% subset(subset = Ensembl %in% gList | Symbol %in% gList, select = c(Cluster, Symbol, Ensembl)) %>%
    mutate(GeneSet = name)
}
rm(name, gList)

for (name in names(Cluster_GOI)) {
  cat("Processing ", name, "\n")
  CB_pHeatmap_methylCluster(name, path = "DataFiles/ClusterHeatmaps/NascentCategories/")
}

df_to_plot <- do.call(rbind, Cluster_GOI) %>% group_by(GeneSet, Cluster) %>% summarise(count = n()) %>%
  mutate(Percent = count / sum(count)) %>% ungroup %>% filter(GeneSet != 'others')
df_to_plot <- Methyl_Clusters %>% group_by(Cluster) %>% summarise(count = n()) %>% mutate(Percent = count / sum(count), GeneSet = 'All') %>% ungroup() %>%
  rbind(df_to_plot, .) %>% mutate(Percent = if_else(Cluster == 'C5', Percent, Percent)) %>%
  mutate(GeneSet = factor(GeneSet, levels = c('upPI_nsGB', 'nsPI_upGB', 'downPI_upGB', 'downPI_nsGB', 'others', 'All')))

#Barplot
ggplot(df_to_plot, aes(x = GeneSet, y = Percent, fill = Cluster)) +
  geom_bar(stat = 'identity', position = position_stack(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0.01, 0), limits = c(0, 1.01))+
  scale_fill_brewer(palette = 'RdYlBu') +
  labs(y = '% Genes', fill = 'Cluster', x = element_blank()) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.8),
        plot.margin = margin(b = 20))

#Piechart
df_to_plot <- subset(df_to_plot, GeneSet != 'All')
ggplot(df_to_plot, aes(x = "", y = Percent, fill = Cluster)) +
  geom_bar(stat = 'identity', width = 1) +
  coord_polar("y") +
  facet_wrap(~GeneSet) +
  scale_fill_brewer(palette = 'RdYlBu') +
  scale_y_continuous(breaks = NULL) +
  labs(x = NULL, y = NULL, fill = 'Clusters') +
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.position = 'right', legend.direction = 'vertical',
        strip.text = element_text(size = 12))


## Plot upregulated genes from full gene DE analysis from PReCIS-seq data
Nascent_full_DEs <- read.xlsx('/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/Precise2024_DiffExpressedGenes_fullGene.xlsx', rowNames = TRUE)

Nascent_up_genes <- subset(Nascent_full_DEs, log2FoldChange > 1 & padj < 0.05, select = 'Ensembl', drop = TRUE)

## Plot nascent gene list with 3 categories, cluster heatmaps
Cluster_GOI <- list()
for (name in names(Nascent_GOI)) {
  gList <- Nascent_GOI[[name]]
  Cluster_GOI[[name]] <- Methyl_Clusters %>% subset(subset = Ensembl %in% gList | Symbol %in% gList, select = c(Cluster, Symbol, Ensembl)) %>%
    mutate(GeneSet = name)
}
rm(name, gList)

Cluster_GOI_3 <- list()
Cluster_GOI_3[['Category_1']] <- Cluster_GOI$upPI_nsGB %>% mutate(GeneSet = 'Category_1')
Cluster_GOI_3[['Category_2']] <- Cluster_GOI$nsPI_upGB %>% mutate(GeneSet = 'Category_2')
Cluster_GOI_3[['Category_3']] <- rbind(Cluster_GOI$downPI_upGB, Cluster_GOI$downPI_nsGB) %>% mutate(GeneSet = 'Category_3')
Cluster_GOI_3[['Ohters']] <- Cluster_GOI$others %>% mutate(GeneSet = 'Others')
Cluster_GOI <- Cluster_GOI_3

for (name in names(Cluster_GOI)) {
  cat("Processing ", name, "\n")
  CB_pHeatmap_methylCluster(name, path = "DataFiles/ClusterHeatmaps/Nascent_3_Categories/")
}

## Plot nascent gene list with 3 categories, promoter region methylation

# Make merged gene set bed files (merging downPI_nsGB and downPI_upGB)
local({
  gr <- subset(all_genes, gene_id %in% Cluster_GOI$Category_3$Ensembl)
  CB_export_bed(gr, fileName = 'DataFiles/Methyl_Plotting/Genebody_Set_downPI_both.bed')
})

local({
  set.seed(6308)
  rdm_df <- sample_n(as.data.frame(all_genes), 150)
  CB_export_bed(rdm_df, 'DataFiles/Methyl_Plotting/GeneBody_random_428.bed', identifier = 'gene_id')
})

command <- paste0("computeMatrix reference-point -S ", "DataFiles/H3K9me3_IgG_subtracted.bw",
                  " -R ", "DataFiles/Methyl_Plotting/Genebody_Set_upPI_nsGB.bed", " ",
                  "DataFiles/Methyl_Plotting/Genebody_Set_nsPI_upGB.bed", " ",
                  "DataFiles/Methyl_Plotting/Genebody_Set_downPI_both.bed", " ",
                  "DataFiles/Methyl_Plotting/GeneBody_random_428.bed",
                  " -o ", "DataFiles/Methyl_Plotting/matrix_nascent_set.gz",
                  " --referencePoint TSS --binSize 5 --sortRegions descend -p 15 -b 500 -a 500")
system(command)
command <- paste("plotHeatmap -m ", "DataFiles/Methyl_Plotting/matrix_nascent_set.gz",
                 " -o ", "DataFiles/Methyl_Plotting/Methylation_aroundTSS_Nascent_3Set.pdf",
                 " --colorMap RdBu_r --samplesLabel H3K9me3 --regionsLabel Category_1 Category_2 Category_3 Random",
                 " --zMin -0.25 --zMax 0.25 --heatmapWidth 8 --heatmapHeight 16 --yMax 0.04")
system(command)



