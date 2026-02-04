library(GenomicFeatures)
library(tidyverse)
library(rtracklayer)
library(openxlsx)

####### Cluster analysis of the whole genome using 50kb flanking region ########
####### Gene body region is scaled to 50kb for plotting purpose, and binSize=100bp for computation
CB_export_bed(all_genes, "DataFiles/ClusterHeatmaps_50kb/Genebody_allGenes.bed", identifier = 'gene_id')

command <- paste0("computeMatrix scale-regions -S ", "DataFiles/H3K9me3_IgG_subtracted.bw",
                  " -R ", "DataFiles/ClusterHeatmaps_50kb/Genebody_allGenes.bed",
                  " -o ", "DataFiles/ClusterHeatmaps_50kb/matrix_all_genes_50kbFlanking.gz",
                  " --regionBodyLength 10000 --binSize 100 --sortRegions descend -p 15 -b 10000 -a 10000 &")
system(command)
system("gunzip DataFiles/ClusterHeatmaps_50kb/matrix_all_genes_50kbFlanking.gz")
methyl_matrix <- read.table("DataFiles/ClusterHeatmaps_50kb/matrix_all_genes_50kbFlanking", skip = 1, sep = '\t')
methyl_matrix <- na.omit(methyl_matrix)
methyl_meta <- methyl_matrix[,1:6]
methyl_data <- methyl_matrix[,7:ncol(methyl_matrix)]

methyl_data <- scale(methyl_data)
kmeans_result <- kmeans(methyl_data, centers = 5)
Methyl_Clusters_50kb <- cbind(methyl_meta, kmeans_result$cluster)
colnames(Methyl_Clusters_50kb)[7] <- 'Cluster'
table(Methyl_Clusters_50kb$Cluster)

names(Methyl_Clusters_50kb) <- c('Chrom', 'ChromStart', 'ChromEnd', 'Ensembl', 'Score', 'Strand', 'Cluster')
gene_name <- info[,2:3] %>% rename(Ensembl = GENEID, Name = GENENAME)
Methyl_Clusters_50kb <- left_join(Methyl_Clusters_50kb, gene_name, by = 'Ensembl')
rm(gene_name)

#Switch clusters names in an easy-to-interpret order [effectively rank gene number for low to high]
ClusterNames <- c("C4", "C1", "C5", "C2", "C3")
Methyl_Clusters_50kb$Cluster <- ClusterNames[Methyl_Clusters_50kb$Cluster]

write.xlsx(Methyl_Clusters_50kb, "DataFiles/ClusterHeatmaps_50kb/Methyl_Clusters_50kb.xlsx")
rm(Methyl_Clusters_50kb)

Methyl_Clusters_50kb <- read.xlsx("DataFiles/ClusterHeatmaps_50kb/Methyl_Clusters_50kb.xlsx", colNames = TRUE)
table(Methyl_Clusters_50kb$Cluster)


rm(methyl_matrix, methyl_meta, methyl_data, kmeans_result)

######## Plotting heatmap for the clusters obtained ######
#subset total of 1000 genes, proportionally distributed among C1-C4 clusters; Subset 200 genes from C5
subsetGenes_clusterred <- subset(Methyl_Clusters_50kb, Cluster != 'C5')
set.seed(428)
subsetGenes_clusterred <- subsetGenes_clusterred %>% sample_n(size = 1000)

for (set in unique(subsetGenes_clusterred$Cluster)) {
  Genes <- subset(subsetGenes_clusterred, Cluster == set, select = Ensembl, drop = TRUE)
  gr <- subset(all_genes, gene_id %in% Genes)
  CB_export_bed(gr, paste0('DataFiles/ClusterHeatmaps_50kb/GB_Genes_', set, '.bed'))
}
rm(set, Genes, gr)

local({
  set.seed(428)
  C5_cluster <- subset(Methyl_Clusters_50kb, Cluster == 'C5') %>% sample_n(size = 200)
  gr <- subset(all_genes, gene_id %in% C5_cluster$Ensembl)
  CB_export_bed(gr, paste0('DataFiles/ClusterHeatmaps_50kb/GB_Genes_C5.bed'))
})

system2("deeptools", args = "--version")
command <- paste0("computeMatrix scale-regions -S ", "DataFiles/H3K9me3_IgG_subtracted.bw",
                  " -R ", "DataFiles/ClusterHeatmaps_50kb/GB_Genes*.bed",
                  " -o ", "DataFiles/ClusterHeatmaps_50kb/matrix_GB_Clusters_SelectedGenes.gz",
                  " --regionBodyLength 50000 --binSize 100 --sortRegions descend -p 15 -b 50000 -a 50000")
system(command)
command <- paste("plotHeatmap -m ", "DataFiles/ClusterHeatmaps_50kb/matrix_GB_Clusters_SelectedGenes.gz",
                 " -o ", "DataFiles/ClusterHeatmaps_50kb/GB_Clusters_SelectedGenes_C1toC5.pdf",
                 " --colorMap RdBu_r --samplesLabel H3K9me3 --regionsLabel C1 C2 C3 C4 C5",
                 " --zMin -0.4 --zMax 1 --yMin -0.3 yMax 2 --heatmapWidth 8 --heatmapHeight 16 --dpi 150")
system(command)

rm(subsetGenes_clusterred)
file.remove('DataFiles/ClusterHeatmaps_50kb/matrix_GB_Clusters_SelectedGenes.gz')
file.remove(list.files("DataFiles/ClusterHeatmaps_50kb/", full.names = TRUE)[grep(".bed", list.files("DataFiles/ClusterHeatmaps/"))])

######## Compare the gene list with that obtained using 10kb flanking region ######
library(ggvenn)
library(patchwork)

Methyl_Clusters_10kb <- read.xlsx("DataFiles/Methyl_Annotation/Methyl_clusters.xlsx", colNames = TRUE)
table(Methyl_Clusters_10kb$Cluster)
table(Methyl_Clusters_50kb$Cluster)

Genes_in_Cluster_10kb <- split(Methyl_Clusters_10kb$Ensembl, Methyl_Clusters_10kb$Cluster)
Genes_in_Cluster_50kb <- split(Methyl_Clusters_50kb$Ensembl, Methyl_Clusters_50kb$Cluster)

plots <- list()
for (cluster in c('C1', 'C2', 'C3', 'C4', 'C5')){
  lt <- list('10kb' = Genes_in_Cluster_10kb[[cluster]],
             '50kb' = Genes_in_Cluster_50kb[[cluster]])
  plots[[cluster]] <- ggvenn(lt, auto_scale = TRUE)
}

rm(Genes_in_Cluster_10kb, Genes_in_Cluster_50kb)

plots[['C1']] + plots[['C2']] + plots[['C3']] + plots[['C4']] + plots[['C5']] +
  plot_layout(ncol = 3)
