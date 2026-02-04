library(tidyverse)
library(openxlsx)



#Sort input files as required by bedtools closest
CB_export_bed(all_genes, "DataFiles/Methyl_Annotation/GB_all_genes.bed", identifier = 'gene_id')
system("sort -k1,1 -k2,2n DataFiles/Methyl_Annotation/GB_all_genes.bed > DataFiles/Methyl_Annotation/GB_all_genes_sorted.bed")
system("sort -k1,1 -k2,2n DataFiles/Epic2_peaks.bed > DataFiles/Epic2_peaks_sorted.bed")

#Find the closest peak for each gene (peak to gene distance)
command <- paste0("bedtools closest -a ", "DataFiles/Methyl_Annotation/GB_all_genes_sorted.bed",
                  " -b ", "DataFiles/Epic2_peaks_sorted.bed",
                  " -d -t first > ",
                  "DataFiles/Methyl_Annotation/ClosestEpic2Peak_allGenes.bed")
system(command)

Peak_distance <- read.table("DataFiles/Methyl_Annotation/ClosestEpic2Peak_allGenes.bed", sep = "\t") %>%
  select(c(V4, V11, V13)) %>% rename(GeneID = V4, PeakScore = V11, Distance = V13)
Peak_distance <- Peak_distance %>% mutate(Group = case_when(GeneID %in% RNAseq_E16_upGenes ~ "Upregulated",
                                                            GeneID %in% random_Genes ~ "Random",
                                                            .default =  "Others"))
Dist_summary <- Peak_distance %>% subset(Group != 'Others') %>% group_by(Group) %>%
  summarise(MeanDist = mean(Distance), Dist_SE = sd(Distance)/sqrt(n())) %>% ungroup()

ggplot(Dist_summary, aes(x = Group, y = MeanDist)) +
  geom_bar(stat = 'identity', fill = 'darkgreen', width = 0.8, color = 'black', size = 1) +
  geom_errorbar(aes(ymin = MeanDist - Dist_SE, ymax = MeanDist + Dist_SE), width = 0.2, size =1) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = element_blank(), y = "Mean distance to closest H3K9me3 peaks") +
  theme(axis.text.x = element_text(size = 12))


#Analyze nearest peak to promoter regions

#Using Annotation from Gopal
Annotations <- read.table("/Users/kb744/Tumbar Lab Local/Precise-seq_2024/ExpressionAnalysis/RefSeq_MostActiveTranscriptPerGene_110123.bed", 
                          header=FALSE, stringsAsFactors = FALSE)
colnames(Annotations) <- c("chrom", "start", "end", "width", "strand", "transcript_id", "gene_name", "gene_type", "gene_id")

all_genes_GC <-  sort(GRanges(Annotations))
rm(Annotations)

TSS_region <- promoters(all_genes, upstream = 500, downstream = 500, use.names = TRUE)
#TSS_region <- promoters(all_genes, upstream = 500, downstream = 500, use.names = TRUE)
TSS_region <- subset(TSS_region, start(TSS_region) >= 0)
CB_export_bed(TSS_region, "DataFiles/Methyl_Annotation/PR_all_genes.bed", identifier = 'gene_id')
system("sort -k1,1 -k2,2n DataFiles/Methyl_Annotation/PR_all_genes.bed > DataFiles/Methyl_Annotation/PR_all_genes_sorted.bed")
command <- paste0("bedtools closest -a ", "DataFiles/Methyl_Annotation/PR_all_genes_sorted.bed",
                  " -b ", "DataFiles/Epic2_peaks_sorted.bed",
                  " -D a -t first > ",
                  "DataFiles/Methyl_Annotation/ClosestEpic2Peak_toPromoters.bed")
system(command)

## Analyze distance with bulk RNAseq data
RNAseq_E16_DEs <- read.xlsx('/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/RNAseq_E16_genotypeOnly_noH614_DEgenes_annotated.xlsx',
                            sheet = 1)
Peak_distance <- read.table("DataFiles/Methyl_Annotation/ClosestEpic2Peak_toPromoters.bed", sep = "\t") %>%
  select(c(V4, V10, V11, V13)) %>% rename(ensembl = V4, PeakScore = V11, Distance = V13, PeakID = V10) %>% filter(PeakScore > 0) %>%
  merge(RNAseq_E16_DEs, ., by = 'ensembl', all = FALSE) %>%
  mutate(Group = if_else(log2FoldChange > 1 & padj < 0.05, "Upregulated", "Others")) %>%
  mutate(Group = factor(Group, levels = c('Upregulated', 'Others'))) %>%
  mutate(LogDistance = case_when(Distance > 0 ~ log10(Distance),
                                 Distance == 0 ~ 0,
                                 Distance < 0 ~ -log10(-Distance))) %>%
  arrange(desc(Group))

ggplot(Peak_distance, aes(x = Distance, y = log2FoldChange, color = Group)) +
  geom_point(size = 1) +
  scale_color_manual(values = c('red', 'gray')) +
  labs(x = 'Distance from TSS to nearest H3K9me3 peak', y = 'RNA Log2(TKO/Ctrl)', color = NULL)

ggplot(Peak_distance, aes(x = Group, y = abs(Distance))) +
  geom_bar(stat = 'summary', fun = 'mean', fill = 'darkgreen', width = 0.8, color = 'black', size = 1) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, size =1) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = element_blank(), y = "Mean distance") +
  theme(axis.text.x = element_text(size = 12))

ggplot(Peak_distance, aes(x = abs(LogDistance), color = Group)) +
  stat_ecdf() +
  scale_color_manual(values = c('red', 'black')) +
  labs(x = 'Log Distance from TSS to nearest H3K9me3 peak', y = 'Cumulative Fraction') +
  theme(axis.title = element_text(size = 10))



#### Analyze nearest peak to promoter for nascent gene categories
PreciseSeqFC_table <- read.xlsx("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/Precise2024_FoldChangeTable_annotated.xlsx")
PreciseSeqFC_table$Category <- factor(PreciseSeqFC_table$Category,
                                      levels = c('upPI_nsGB', 'nsPI_upGB', 'downPI_upGB', 'downPI_nsGB', 'others'))

Peak_distance <- read.table("DataFiles/Methyl_Annotation/ClosestEpic2Peak_toPromoters.bed", sep = "\t") %>%
  select(c(V4, V11, V13)) %>% rename(GeneID = V4, PeakScore = V11, Distance = V13) %>%
  filter(PeakScore > 0) %>%
  inner_join(PreciseSeqFC_table[,c('GeneID', 'Symbol', 'GENEBIOTYPE', 'DEseq_log2FC_PI', 'DEseq_log2FC_GB', 'Category')], by = 'GeneID') %>%
  arrange(desc(Category))

df_to_plot <- Peak_distance
ggplot(df_to_plot, aes(x = Distance, y = DEseq_log2FC_PI, color = Category)) +
  geom_point(size = 1) +
  scale_color_manual(values = c('red', 'orange', 'lightblue', 'yellowgreen', 'gray')) +
  labs(x = 'Distance from TSS to nearest H3K9me3 peak', y = 'PI log2(TKO/Ctrl)', color = NULL)
ggplot(df_to_plot, aes(x = Distance, y = DEseq_log2FC_GB, color = Category)) +
  geom_point(size = 1) +
  scale_color_manual(values = c('red', 'orange', 'lightblue', 'yellowgreen', 'gray')) +
  labs(x = 'Distance from TSS to nearest H3K9me3 peak', y = 'GB log2(TKO/Ctrl)', color = NULL)

ggplot(df_to_plot, aes(x = Category, y = abs(Distance))) +
  geom_bar(stat = 'summary', fun = 'mean', fill = 'darkgreen', width = 0.8, color = 'black', size = 1) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, size =1) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = element_blank(), y = "Mean distance of closest H3K9me3 peaks to TSS") +
  theme(axis.text.x = element_text(size = 12))

df_summary <- df_to_plot %>% group_by(Category) %>% summarise(Pct_Internal = mean(Distance == 0 ))
ggplot(df_summary, aes(x = Category, y = Pct_Internal)) +
  geom_bar(stat = 'identity', fill = 'darkgreen', width = 0.8, color = 'black', size = 1) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(scale = 100)) +
  labs(x = element_blank(), y = 'Percent of genes with TSS methylation')

df_for_ecdf <- df_to_plot %>% subset(Category != 'others') %>%
  mutate(LogDistance = case_when(Distance < 0 ~ -log10(-Distance),
                                 Distance == 0 ~ 0,
                                 Distance > 0 ~ log10(Distance)))
ggplot(df_for_ecdf, aes(x = abs(LogDistance), color = Category)) +
  stat_ecdf()


###Distance analysis for scRNAseq related genesets

#Get gene lists related to scRNAseq
local({
  wb <- loadWorkbook("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/scRNA_DE_genes_v4_annotated.xlsx")
  ClusterDEgenes <<- list()
  for (name in names(wb)) {
    ClusterDEgenes[[name]] <<- readWorkbook(wb, sheet = name)
  }
})

Cluster_idents <- c('0' = 'Basal_1', '1' = 'Diff_1', '2' = 'Basal_2', '3' = 'HairFollicle',
                    '4' = 'AberrBasal', '5' = 'Basal_3', '6' = 'Diff_2', '7' = 'Diff_3',
                    '8' = 'Merkel', '9' = 'AberrDiff')
names(ClusterDEgenes) <- Cluster_idents

Peak_distance <- read.table("DataFiles/Methyl_Annotation/ClosestEpic2Peak_toPromoters.bed", sep = "\t") %>%
  select(c(V4, V11, V13)) %>% rename(GeneID = V4, PeakScore = V11, Distance = V13) %>% column_to_rownames('GeneID')

scRNA_upGene <- lapply(ClusterDEgenes[-9], function(dt){
  upTable <- subset(dt, p_val_adj < 0.05, abs(avg_log2FC) > 1) %>%
    mutate(GENEID = if_else(is.na(GENEID), GeneID, GENEID))
  Pd <- Peak_distance[upTable$GENEID,]
  upTable$Distance <- Pd$Distance
  upTable$Score <- Pd$PeakScore
  return(upTable)
})

df_to_plot <- do.call(rbind, scRNA_upGene) %>% na.omit()
df_to_plot$Cluster <- sub("\\..*", "", rownames(df_to_plot))
rownames(df_to_plot) <- NULL
df_to_plot <- df_to_plot %>%
  mutate(Group = case_when(Cluster %in% c('AberrBasal', 'AberrDiff', 'HairFollicle') ~ Cluster,
                           Cluster %in% c('Basal_1', 'Basal_2', 'PrePlacode') ~ 'Basal',
                           Cluster %in% c('Diff_1', 'Diff_2', 'Diff_3') ~ 'Diff')) %>%
  filter(Group != 'HairFollicle')

library(scales)
ggplot(df_to_plot, aes(x = Distance, y = avg_log2FC)) +
  geom_point() +
  scale_x_continuous(labels = label_number(scale_cut = cut_si(''))) +
  facet_wrap(~Group, nrow = 1) +
  labs(x = 'Distance from TSS to nearest H3K9me3 peak',
       y = 'Log2(TKO/Ctrl)')

ggplot(df_to_plot, aes(x = Group, y = abs(Distance))) +
  geom_bar(stat = 'summary', fun = 'mean', fill = 'darkgreen', width = 0.8, color = 'black', size = 1) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, size =1) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = element_blank(), y = "Mean distance to closest H3K9me3 peaks") +
  theme(axis.text.x = element_text(size = 12))



