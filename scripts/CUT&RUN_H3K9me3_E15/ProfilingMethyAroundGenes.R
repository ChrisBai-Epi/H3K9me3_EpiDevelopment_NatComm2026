library(GenomicFeatures)
library(tidyverse)
library(rtracklayer)
library(openxlsx)


#Get and define gene region, promoter region and flanking regions
txdb_ensembl <- makeTxDbFromEnsembl(organism = 'Mus musculus', release = 110)
all_genes <- genes(txdb_ensembl)
info <- select(Ens_mm39_110, keys = unlist(all_genes$gene_id), keytype = 'GENEID', 
               columns = c('GENEBIOTYPE', 'GENEID', 'GENENAME', 'SYMBOL', 'DESCRIPTION'))
all_genes$Biotype <- info$GENEBIOTYPE
all_genes$Symbol <- info$SYMBOL
rm(info)

chrs <- c(paste0('chr', 1:19), 'chrX', 'chrY', 'chrM')
seqlevelsStyle(all_genes) <- 'UCSC'
all_genes <- keepSeqlevels(all_genes, chrs, pruning.mode = 'tidy')
all_genes <- sort(all_genes)
names(all_genes) <- NULL
all_promoters <- promoters(all_genes, upstream = 500, downstream = 500, use.names = TRUE)
all_upstream_10k <- flank(all_genes, width = 10000, start = TRUE, use.names = TRUE)
all_downstream_10k <- flank(all_genes, width = 10000, start = FALSE, use.names = TRUE)


#######Export regions as bed files
CB_export_bed <- function(gr, fileName, identifier = 'Symbol') {
  df <- as.data.frame(gr)
  df$start <- df$start - 1
  df$Score <- 0
  bed_df <- df[, c('seqnames', 'start', 'end', identifier, 'Score', 'strand', 'gene_id', 'Symbol', 'Biotype')]
  write.table(bed_df, file = fileName, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

CB_export_bed(all_genes, 'DataFiles/all_genes.bed')
CB_export_bed(all_promoters, 'DataFiles/all_promoters.bed')
CB_export_bed(all_upstream_10k, 'DataFiles/all_upstream_10k.bed')
CB_export_bed(all_downstream_10k, 'DataFiles/all_downstream_10k.bed')

gr_up_DEgenes <- all_genes[all_genes$gene_id %in% Up_genes_combined | all_genes$Symbol %in% Up_genes_combined]
gr_down_DEgenes<- all_genes[all_genes$gene_id %in% Down_genes_combined | all_genes$Symbol %in% Down_genes_combined]
gr_up_TKO_markers <- all_genes[all_genes$gene_id %in% TKO_Markers_Up | all_genes$Symbol %in% TKO_Markers_Up]
gr_down_TKO_markers <- all_genes[all_genes$gene_id %in% TKO_Markers_Down | all_genes$Symbol %in% TKO_Markers_Down]

CB_export_bed(gr_up_DEgenes, 'DataFiles/Regions_of_interest/gr_up_DEgenes.bed')
CB_export_bed(gr_down_DEgenes, 'DataFiles/Regions_of_interest/gr_down_DEgenes.bed')
CB_export_bed(gr_up_TKO_markers, 'DataFiles/Regions_of_interest/gr_up_TKO_markers.bed')
CB_export_bed(gr_down_TKO_markers, 'DataFiles/Regions_of_interest/gr_down_TKO_markers.bed')



######Plotting Heatmaps of methylation around genes

#10kb around TSS
all_promoters_10k <- promoters(all_genes, upstream = 10000, downstream = 10000, use.names = TRUE)
all_promoters_10k <- GenomicRanges::sort(all_promoters_10k, ignore.strand = TRUE)
all_genes <- GenomicRanges::sort(all_genes, ignore.strand = TRUE)


#### Plot with deeptools
#### Import cluster information obtained from deeptools

methyl_cluster_scUp <- read.table("DataFiles/Methyl_Annotation/GB_scUp_H3K9me3_kmeans_5.bed", sep = '\t')
methyl_cluster_bulk <- read.table("DataFiles/Methyl_Annotation/GB_RNAseq_Up_H3K9me3_kmeans_5.bed", sep = '\t')
methyl_cluster_nascent <- read.table("DataFiles/Methyl_Annotation/GB_Nascent_Up_H3K9me3_kmeans_5.bed", sep = '\t')

methyl_cluster_all <- methyl_cluster_scUp[, c(4, 13)] %>%
  merge(x=., y = methyl_cluster_bulk[,c(4,13)], by = 'V4', all = TRUE) %>%
  merge(x=., y = methyl_cluster_nascent[,c(4,13)], by = 'V4', all = TRUE)
colnames(methyl_cluster_all) <- c('GeneSymbol', 'scRNAseq', 'BulkRNA', 'NascentRNA')
rm(methyl_cluster_bulk, methyl_cluster_scUp, methyl_cluster_nascent)

write.xlsx(methyl_cluster_all, "/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/MethylationClusterInfo.xlsx")
addGeneInfo("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/MethylationClusterInfo.xlsx", enbl)
#This part (above) might not be super useful, since we got manual kmeans analysis below

#######Cluster analysis of the whole genome
CB_export_bed(all_genes, "DataFiles/ClusterHeatmaps/Genebody_allGenes.bed", identifier = 'gene_id')

command <- paste0("computeMatrix scale-regions -S H3K9me3_IgG_subtracted.bw -R DataFiles/ClusterHeatmaps/Genebody_allGenes.bed ",
                  "-o ClusterHeatmaps/matrix_GB_allGenes.gz --regionBodyLength 10000 --binSize 20 --sortRegions descend -p 15 -b 10000 -a 10000 &")
system(command)
system("gunzip ClusterHeatmaps/matrix_GB_allGenes.gz")
methyl_matrix <- read.table("ClusterHeatmaps/matrix_GB_allGenes", skip = 1, sep = '\t')
methyl_matrix <- na.omit(methyl_matrix)
methyl_meta <- methyl_matrix[,1:6]
methyl_data <- methyl_matrix[,7:ncol(methyl_matrix)]

methyl_data <- scale(methyl_data)
kmeans_result <- kmeans(methyl_data, centers = 5)
Methyl_Clusters <- cbind(methyl_meta, kmeans_result$cluster)
colnames(Methyl_Clusters)[7] <- 'Cluster'
table(Methyl_Clusters$Cluster)

names(Methyl_Clusters) <- c('Chrom', 'ChromStart', 'ChromEnd', 'Ensembl', 'Score', 'Strand', 'Cluster')
Methyl_Clusters$Symbol <- info[Methyl_Clusters$Ensembl,]$SYMBOL

#Switch clusters names in an easy-to-interpret order
ClusterNames <- c("C1", "C5", "C2", "C3", "C4")
Methyl_Clusters$Cluster <- ClusterNames[Methyl_Clusters$Cluster]
rm(Methyl_Clusters)

write.xlsx(Methyl_Clusters, "DataFiles/Methyl_Annotation/Methyl_clusters.xlsx")
Methyl_Clusters <- read.xlsx("DataFiles/Methyl_Annotation/Methyl_clusters.xlsx", colNames = TRUE)
table(Methyl_Clusters$Cluster)

rm(methyl_matrix, methyl_meta, methyl_data, kmeans_result)

#Subset genes for heatmap plotting
#subset equal number of genes from each cluster
set.seed(428)
subsetGenes_clusterred <- Methyl_Clusters %>% group_by(Cluster) %>%
  sample_n(size = min(200, n())) %>% ungroup()

#subset total of 1000 genes, proportionally distributed among the clusters
set.seed(428)
subsetGenes_clusterred <- Methyl_Clusters %>% sample_n(size = 1000)

for (set in unique(subsetGenes_clusterred$Cluster)) {
  Genes <- subset(subsetGenes_clusterred, Cluster == set, select = Ensembl, drop = TRUE)
  gr <- subset(all_genes, gene_id %in% Genes)
  CB_export_bed(gr, paste0('DataFiles/ClusterHeatmaps/GB_Genes_', set, '.bed'))
}
rm(set, Genes, gr)

system2("deeptools", args = "--version")
command <- paste0("computeMatrix scale-regions -S ", "DataFiles/H3K9me3_IgG_subtracted.bw",
                  " -R ", "DataFiles/ClusterHeatmaps/GB_Genes*.bed",
                  " -o ", "DataFiles/ClusterHeatmaps/matrix_GB_Clusters_SelectedGenes.gz",
                  " --regionBodyLength 10000 --binSize 20 --sortRegions descend -p 15 -b 10000 -a 10000")
system(command)
command <- paste("plotHeatmap -m ", "DataFiles/ClusterHeatmaps/matrix_GB_Clusters_SelectedGenes.gz",
                 " -o ", "DataFiles/ClusterHeatmaps/GB_Clusters_SelectedGenes_200perCluster.pdf",
                 " --colorMap RdBu_r --samplesLabel H3K9me3 --regionsLabel C1 C2 C3 C4 C5 ",
                 "--zMin -0.4 --zMax 1 --yMin -0.3 yMax 2 --heatmapWidth 8 --heatmapHeight 16")
system(command)

rm(subsetGenes_clusterred)
file.remove('DataFiles/ClusterHeatmaps/matrix_GB_Clusters_SelectedGenes.gz')
file.remove(list.files("DataFiles/ClusterHeatmaps/", full.names = TRUE)[grep(".bed", list.files("DataFiles/ClusterHeatmaps/"))])

#subset total of 1000 genes, proportionally distributed among C1-C4 clusters; Subset 200 genes from C5
subsetGenes_clusterred <- subset(Methyl_Clusters, Cluster != 'C5')
set.seed(428)
subsetGenes_clusterred <- subsetGenes_clusterred %>% sample_n(size = 1000)

for (set in unique(subsetGenes_clusterred$Cluster)) {
  Genes <- subset(subsetGenes_clusterred, Cluster == set, select = Ensembl, drop = TRUE)
  gr <- subset(all_genes, gene_id %in% Genes)
  CB_export_bed(gr, paste0('DataFiles/ClusterHeatmaps/GB_Genes_', set, '.bed'))
}
rm(set, Genes, gr)

local({
  set.seed(428)
  C5_cluster <- subset(Methyl_Clusters, Cluster == 'C5') %>% sample_n(size = 200)
  gr <- subset(all_genes, gene_id %in% C5_cluster$Ensembl)
  CB_export_bed(gr, paste0('DataFiles/ClusterHeatmaps/GB_Genes_C5.bed'))
})

system2("deeptools", args = "--version")
command <- paste0("computeMatrix scale-regions -S ", "DataFiles/H3K9me3_IgG_subtracted.bw",
                  " -R ", "DataFiles/ClusterHeatmaps/GB_Genes*.bed",
                  " -o ", "DataFiles/ClusterHeatmaps/matrix_GB_Clusters_SelectedGenes.gz",
                  " --regionBodyLength 10000 --binSize 20 --sortRegions descend -p 15 -b 10000 -a 10000")
system(command)
command <- paste("plotHeatmap -m ", "DataFiles/ClusterHeatmaps/matrix_GB_Clusters_SelectedGenes.gz",
                 " -o ", "DataFiles/ClusterHeatmaps/GB_Clusters_SelectedGenes_C1toC5.pdf",
                 " --colorMap RdBu_r --samplesLabel H3K9me3 --regionsLabel C1 C2 C3 C4 C5",
                 " --zMin -0.4 --zMax 1 --yMin -0.3 yMax 2 --heatmapWidth 8 --heatmapHeight 16 --dpi 150")
system(command)

rm(subsetGenes_clusterred)
file.remove('DataFiles/ClusterHeatmaps/matrix_GB_Clusters_SelectedGenes.gz')
file.remove(list.files("DataFiles/ClusterHeatmaps/", full.names = TRUE)[grep(".bed", list.files("DataFiles/ClusterHeatmaps/"))])



###Profile distribution of clusters across different gene sets
##Getting different gene sets
RNAseq_E15_DEs <- read.xlsx('/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/RNAseq_E15_TKOvsCtrl_annotated.xlsx')
RNAseq_E16_DEs <- read.xlsx('/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/RNAseq_E16_genotypeOnly_noH614_DEgenes_annotated.xlsx',
                            sheet = 2)

RNAseq_E15_upGenes <- subset(RNAseq_E15_DEs, log2FoldChange > 0 & padj < 0.05, select = 'ensembl', drop = TRUE)
RNAseq_E16_upGenes <- subset(RNAseq_E16_DEs, log2FoldChange > 0, select = 'ensembl', drop = TRUE)
RNAseq_E16_downGenes <- subset(RNAseq_E16_DEs, log2FoldChange < 0, select = 'ensembl', drop = TRUE)

set.seed(428)
random_Genes <- sample_n(as.data.frame(all_genes), 300)
random_Genes <- random_Genes$gene_id

list_GOI <- list('SC_upGenes' = Up_genes_combined, 'SC_downGenes' = Down_genes_combined, 'Nascent_up' = Nascent_up_genes,
                 'RNAseq_up' = RNAseq_E16_upGenes, 'RNAseq_down' = RNAseq_E16_downGenes, 'Random' = random_Genes)

#Methyl_Clusters$Cluster <- paste0("C", Methyl_Clusters$Cluster)
Cluster_GOI <- list()
for (name in names(list_GOI)) {
  gList <- list_GOI[[name]]
  Cluster_GOI[[name]] <- Methyl_Clusters %>% subset(subset = Ensembl %in% gList | Symbol %in% gList, select = c(Cluster, Symbol, Ensembl)) %>%
    mutate(GeneSet = name)
}
rm(name, gList)

df_to_plot <- do.call(rbind, Cluster_GOI) %>% group_by(GeneSet, Cluster) %>% summarise(count = n()) %>%
  mutate(Percent = count / sum(count)) %>% ungroup
df_to_plot <- Methyl_Clusters %>% group_by(Cluster) %>% summarise(count = n()) %>%
  mutate(Percent = count / sum(count), GeneSet = 'Genome') %>% ungroup() %>%
  rbind(df_to_plot, .)
df_to_plot$GeneSet <- factor(df_to_plot$GeneSet, levels = c("Genome", "RNAseq_up", "RNAseq_down", 'Nascent_up', "SC_upGenes", "SC_downGenes", "Random"))
df_to_plot <- subset(df_to_plot, GeneSet %in% c('Genome', 'RNAseq_up'))

# Barplots
sub_df_to_plot <- df_to_plot %>% subset(GeneSet == 'Genome')
p1 <- ggplot(sub_df_to_plot, aes(x = GeneSet, y = Percent, fill = Cluster)) +
  geom_bar(stat = 'identity', position = position_stack(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0.005))+
  scale_fill_brewer(palette = 'RdYlBu') +
  labs(y = '% Genes', fill = 'Cluster', x = element_blank()) +
  theme(legend.position = 'none')

sub_df_to_plot <- df_to_plot %>% subset(GeneSet == 'RNAseq_up')
p2 <- ggplot(sub_df_to_plot, aes(x = GeneSet, y = Percent, fill = Cluster)) +
  geom_bar(stat = 'identity', position = position_stack(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0.005))+
  scale_fill_brewer(palette = 'RdYlBu') +
  labs(y = '% Genes', fill = 'Cluster', x = element_blank())

p1 | p2

rm(sub_df_to_plot)

# Piecharts
ggplot(df_to_plot, aes(x = "", y = Percent, fill = Cluster)) +
  geom_bar(stat = 'identity', width = 1) +
  coord_polar("y") +
  facet_wrap(~GeneSet) +
  scale_fill_brewer(palette = 'RdYlBu') +
  scale_y_continuous(breaks = NULL) +
  labs(x = NULL, y = NULL, fill = 'Clusters') +
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.position = 'bottom', legend.direction = 'horizontal',
        strip.text = element_text(size = 12))

## Plot heatmap around promoter regions
# For Nascent_up genes
local({
  gr <- subset(all_genes, gene_id %in% Nascent_up_genes)
  CB_export_bed(gr, 'DataFiles/Methyl_Plotting/GeneBody_Nascent_up.bed', identifier = 'gene_id')
})

local({
  set.seed(9854)
  rdm_df <- sample_n(as.data.frame(all_genes), 300)
  CB_export_bed(rdm_df, 'DataFiles/Methyl_Plotting/GeneBody_random_testing.bed', identifier = 'gene_id')
})

command <- paste0("computeMatrix reference-point -S ", "DataFiles/H3K9me3_IgG_subtracted.bw",
                  " -R ", "DataFiles/Methyl_Plotting/GeneBody_Nascent_up.bed", " ", "DataFiles/Methyl_Plotting/GeneBody_random_testing.bed",
                  " -o ", "DataFiles/Methyl_Plotting/matrix_nascent_up.gz",
                  " --referencePoint TSS --binSize 100 --sortRegions descend -p 15 -b 10000 -a 10000")
system(command)
command <- paste("plotHeatmap -m ", "DataFiles/Methyl_Plotting/matrix_nascent_up.gz",
                 " -o ", "DataFiles/Methyl_Plotting/Methylation_aroundTSS_Nascent_Up.pdf",
                 " --colorMap RdBu_r --samplesLabel H3K9me3 --regionsLabel Nascent_up Random",
                 " --zMin -0.25 --zMax 0.25 --heatmapWidth 8 --heatmapHeight 16")
system(command)

#For bulk_up genes
local({
  gr <- subset(all_genes, gene_id %in% RNAseq_E16_upGenes)
  CB_export_bed(gr, 'DataFiles/Methyl_Plotting/GeneBody_RNAseq_up.bed', identifier = 'gene_id')
})


command <- paste0("computeMatrix reference-point -S ", "DataFiles/H3K9me3_IgG_subtracted.bw",
                  " -R ", "DataFiles/Methyl_Plotting/GeneBody_RNAseq_up.bed", " ", "DataFiles/Methyl_Plotting/Genebody_Random.bed",
                  " -o ", "DataFiles/Methyl_Plotting/matrix_RNAseq_up.gz",
                  " --referencePoint TSS --binSize 100 --sortRegions descend -p 15 -b 10000 -a 10000")
system(command)
command <- paste("plotHeatmap -m ", "DataFiles/Methyl_Plotting/matrix_RNAseq_up.gz",
                 " -o ", "DataFiles/Methyl_Plotting/Methylation_aroundTSS_RNAseq_Up.pdf",
                 " --colorMap RdBu_r --samplesLabel H3K9me3 --regionsLabel RNAseq_up Random",
                 " --zMin -0.25 --zMax 0.25 --heatmapWidth 8 --heatmapHeight 16")
system(command)

#Plot methylation around TSS region for gene categories define with PReCIS-seq data
local({
  set.seed(6308)
  rdm_df <- sample_n(as.data.frame(all_genes), 150)
  CB_export_bed(rdm_df, 'DataFiles/Methyl_Plotting/GeneBody_random_428.bed', identifier = 'gene_id')
})

Peak_distance <- read.table("DataFiles/Methyl_Annotation/ClosestEpic2Peak_allGenes.bed", sep = "\t") %>%
  select(c(V4, V11, V13)) %>% rename(GeneID = V4, PeakScore = V11, Distance = V13)
methylated_TSS <- subset(Peak_distance, Distance == 0 & PeakScore > 300)
set.seed(6308)
methylated_TSS <- sample_n(methylated_TSS, 150)
methylated_TSS <- subset(all_genes, gene_id %in% methylated_TSS$GeneID)
CB_export_bed(methylated_TSS, 'DataFiles/Methyl_Plotting/GeneBody_TSSmethylation', identifier = 'gene_id')

command <- paste0("computeMatrix reference-point -S ", "DataFiles/H3K9me3_IgG_subtracted.bw",
                  " -R ", "DataFiles/Methyl_Plotting/Genebody_Set_upPI_nsGB.bed", " ",
                  "DataFiles/Methyl_Plotting/Genebody_Set_nsPI_upGB.bed", " ",
                  "DataFiles/Methyl_Plotting/Genebody_Set_downPI_upGB.bed", " ",
                  "DataFiles/Methyl_Plotting/Genebody_Set_downPI_nsGB.bed", " ",
                  "DataFiles/Methyl_Plotting/GeneBody_random_428.bed", " ",
                  "DataFiles/Methyl_Plotting/GeneBody_TSSmethylation",
                  " -o ", "DataFiles/Methyl_Plotting/matrix_nascent_set.gz",
                  " --referencePoint TSS --binSize 5 --sortRegions descend -p 15 -b 500 -a 500")
system(command)
command <- paste("plotHeatmap -m ", "DataFiles/Methyl_Plotting/matrix_nascent_set.gz",
                 " -o ", "DataFiles/Methyl_Plotting/Methylation_aroundTSS_Nascent_Set.pdf",
                 " --colorMap RdBu_r --samplesLabel H3K9me3 --regionsLabel upPI_nsGB nsPI_upGB downPI_upGB downPI_nsGB random methylated_TSS",
                 " --zMin -0.25 --zMax 0.25 --heatmapWidth 8 --heatmapHeight 16 --yMax 0.08")
system(command)

### Plot Precise-seq upregulated vs unexpresssed genes
Unexpressed <- read_rds('DataFiles/Unexpressed_GeneList.rds')
set.seed(6843)
Unexpressed <- sample(Unexpressed, 250)
Unexpressed <- subset(all_genes, gene_id %in% Unexpressed)
CB_export_bed(Unexpressed, 'DataFiles/Methyl_Plotting/GeneBody_UnexpressedGenes.bed', identifier = 'gene_id')

set.seed(6843)
Unchanged <- sample(Nascent_GOI$others, 250)
Unchanged <- subset(all_genes, gene_id %in% Unchanged)
CB_export_bed(Unexpressed, 'DataFiles/Methyl_Plotting/GeneBody_UnchanegdGenes.bed', identifier = 'gene_id')

rm(Unexpressed, Unchanged)

command <- paste0("computeMatrix reference-point -S ", "DataFiles/H3K9me3_IgG_subtracted.bw",
                  " -R ", "DataFiles/Methyl_Plotting/GeneBody_Nascent_up.bed", " ", "DataFiles/Methyl_Plotting/GeneBody_UnchanegdGenes.bed",
                  " -o ", "DataFiles/Methyl_Plotting/matrix_nascent_up.gz",
                  " --referencePoint TSS --binSize 100 --sortRegions descend -p 15 -b 10000 -a 10000")
system(command)
command <- paste("plotHeatmap -m ", "DataFiles/Methyl_Plotting/matrix_nascent_up.gz",
                 " -o ", "DataFiles/Methyl_Plotting/Methylation_aroundTSS_Nascent_Up_2.pdf",
                 " --colorMap RdBu_r --samplesLabel H3K9me3 --regionsLabel Nascent_up Unchanged",
                 " --zMin -0.25 --zMax 0.25 --heatmapWidth 8 --heatmapHeight 16")
system(command)

#Plot heatmaps for each gene category
for (name in names(Cluster_GOI)) {
  cat("Processing ", name, "\n")
  CB_pHeatmap_methylCluster(name, path = "DataFiles/ClusterHeatmaps/GeneGroups/")
}

###Define the function used for plotting heatmap
CB_pHeatmap_methylCluster <- function(set, path = "./"){
  if (!dir.exists(path)) stop("The provided path is not a valid directory.")
  regionLabel <- character()
  for (cluster in c('C1', 'C2', 'C3', 'C4', 'C5')){
    ids <- subset(Cluster_GOI[[set]], subset = Cluster == cluster, select = Ensembl, drop = TRUE)
    if (length(ids) <= 1) {
      next
    }
    regionLabel <- paste(regionLabel, cluster)
    gr <- subset(all_genes, gene_id %in% ids)
    CB_export_bed(gr, paste0(path, 'Genebody_Set_', cluster, '.bed'), identifier = 'gene_id')
  }
  command <- paste0("computeMatrix scale-regions -S ", "DataFiles/H3K9me3_IgG_subtracted.bw", 
                    " -R ", path, "Genebody_Set* -o ", path, "Matrix_GB_Cluster_Genes.gz", 
                    " --regionBodyLength 10000 --binSize 20 --sortRegions descend -p 15 -b 10000 -a 10000")
  system(command)
  command <- paste0("plotHeatmap -m ", path, "Matrix_GB_Cluster_Genes.gz", 
                    " --colorMap RdBu_r --samplesLabel H3K9me3 --regionsLabel", regionLabel,
                    " --heatmapWidth 8 --heatmapHeight 16 --zMin -0.4 --zMax 1.0 --yMin -0.3 --yMax 2 ", 
                    "-o ", path, "GB_Heatmap_", set, "_Clusters.pdf")
  cat(regionLabel, '\n')
  system(command)
  file.remove(paste0(path, 'Matrix_GB_Cluster_Genes.gz'))
  file.remove(list.files(path, full.names = TRUE)[grep(".bed", list.files(path))])
}

rm(list = ls()[grep("^gr", ls())])
setwd("~/Tumbar Lab Local/CUTandRUN/Functional_Analysis")
save.image()

