library(tidyverse)
library(openxlsx)

##### This script aims to analyze the methylation states for genes with different pausing levels.

### Read in files
files <- list.files(path = "/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/Pausing_Gopal", full.names = TRUE)

GenePausingInfo <- lapply(files, read.xlsx)
GenePausingInfo <- bind_rows(GenePausingInfo)

rm(files)

### Assign clusters
clusterInfo <- Methyl_Clusters[!duplicated(Methyl_Clusters$Symbol), c(8,7)]
rownames(clusterInfo) <- clusterInfo$Symbol

GenePausingInfo$MethylCluster <- clusterInfo[GenePausingInfo$Gene, 2]
df_to_plot <- GenePausingInfo %>% drop_na() %>% group_by(Quadrant, MethylCluster) %>%
  summarise(Count = n()) %>% mutate(Pct = Count / sum(Count)) %>% ungroup()

### Plotting
ggplot(df_to_plot, aes(x = Quadrant, y = Pct, fill = MethylCluster)) +
  geom_bar(stat = 'identity') +
  scale_fill_brewer(palette = 'RdYlBu')

### Plot heatmaps
dataList <- lapply(list.files(path = "/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/Pausing_Gopal", full.names = TRUE),
                   read.xlsx)
Pausing_GOI <- lapply(dataList, function(df)df$Gene)
names(Pausing_GOI) <- c('Expressed_highPaused', 'Expressed_lowPaused', 'Expressed_midPaused',
                        'Unexpressed_highPaused', 'Unexpressed_lowPaused', 'Unexpressed_midPaused')

Cluster_GOI <- list()
for (name in names(Pausing_GOI)) {
  gList <- Pausing_GOI[[name]]
  Cluster_GOI[[name]] <- Methyl_Clusters %>% subset(subset = Ensembl %in% gList | Symbol %in% gList, select = c(Cluster, Symbol, Ensembl)) %>%
    mutate(CellCluster = name)
}
rm(name, gList, dataList)

for (name in names(Cluster_GOI)) {
  cat("Processing ", name, "\n")
  CB_pHeatmap_methylCluster(name, path = "DataFiles/ClusterHeatmaps/PausingCategories/")
}

### Profile promoter region methylation
for (name in names(Pausing_GOI)){
  gr <- subset(all_genes, Symbol %in% Pausing_GOI[[name]])
  CB_export_bed(gr, paste0("DataFiles/Methyl_Plotting/Genebody_Pause_", name, ".bed"), identifier = 'gene_id')
}
rm(gr, name)

local({
  set.seed(6308)
  rdm_df <- sample_n(as.data.frame(all_genes), 1000)
  CB_export_bed(rdm_df, 'DataFiles/Methyl_Plotting/GeneBody_random_1000.bed', identifier = 'gene_id')
})

local({
  Peak_distance <- read.table("DataFiles/Methyl_Annotation/ClosestEpic2Peak_allGenes.bed", sep = "\t") %>%
    select(c(V4, V11, V13)) %>% rename(GeneID = V4, PeakScore = V11, Distance = V13)
  methylated_TSS <- subset(Peak_distance, Distance == 0 & PeakScore > 300)
  set.seed(6308)
  methylated_TSS <- subset(all_genes, gene_id %in% methylated_TSS$GeneID)
  CB_export_bed(methylated_TSS, 'DataFiles/Methyl_Plotting/GeneBody_TSSmethylation_1000', identifier = 'gene_id')
})


command <- paste0("computeMatrix reference-point -S ", "DataFiles/H3K9me3_IgG_subtracted.bw",
                  " -R ", "DataFiles/Methyl_Plotting/Genebody_Pause_Expressed_highPaused.bed", " ",
                  "DataFiles/Methyl_Plotting/Genebody_Pause_Expressed_lowPaused.bed", " ",
                  "DataFiles/Methyl_Plotting/GeneBody_random_1000.bed", " ",
                  "DataFiles/Methyl_Plotting/GeneBody_TSSmethylation_1000", " ",
                  " -o ", "DataFiles/Methyl_Plotting/matrix_Pause_set.gz",
                  " --referencePoint TSS --binSize 5 --sortRegions descend -p 15 -b 500 -a 500")
system(command)
command <- paste("plotHeatmap -m ", "DataFiles/Methyl_Plotting/matrix_Pause_set.gz",
                 " -o ", "DataFiles/Methyl_Plotting/Methylation_aroundTSS_PausingGroups.pdf",
                 " --colorMap RdBu_r --samplesLabel H3K9me3 --regionsLabel Exp_highPause Exp_lowPause Random Methylated",
                 " --zMin -0.25 --zMax 0.25 --heatmapWidth 8 --heatmapHeight 16 --yMax 0.08")
system(command)



