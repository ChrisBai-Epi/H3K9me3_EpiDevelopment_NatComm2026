library(GenomicRanges)
library(tidyverse)
library(patchwork)
library(openxlsx)

####### PART I: Characterizing different types of peaks #######
#### Import Epic2 peak file
Epic2_peaks <- read.table("DataFiles/epic2_peaks_H3K9me3_E15.bed", header = TRUE)
Epic2_peaks$PeakID <- paste0("Epic", row.names(Epic2_peaks))
Epic2_peaks <- Epic2_peaks[, c(1, 2, 3, 11, 5, 6, 7, 8, 9, 10, 4)]
Epic2_peaks <- subset(Epic2_peaks, Chromosome %in% chrs)

write.table(Epic2_peaks[, 1:6], "DataFiles/Epic2_peaks.bed", col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
Epic2_peaks <- GRanges(Epic2_peaks)

summary(width(Epic2_peaks))

#### Compute mean methylation signal over each peak
cmd <- paste0("sort -k1,1 -k2,2n ", "DataFiles/Epic2_peaks.bed", " > ", "DataFiles/Epic2_peaks_sorted.bed")
system(cmd)
cmd <- paste0("sort -k1,1 -k2,2n ", "DataFiles/H3K9me3_IgG_subtracted.bedgraph", " > ", "DataFiles/H3K9me3_IgG_subtracted_sorted.bedgraph")
system(cmd)

cmd <- paste0("bedtools map -c 4 -o mean -a ",
              "DataFiles/Epic2_peaks_sorted.bed",
              " -b ", "DataFiles/H3K9me3_IgG_subtracted.bedgraph",
              " > ", "DataFiles/Epic2_peaks_methylSignalCounted.bed")
system(cmd)

MethylSignal <- read.table("DataFiles/Epic2_peaks_methylSignalCounted.bed") %>% column_to_rownames("V4")
Epic2_peaks$MeanIntensity <- MethylSignal[Epic2_peaks$PeakID, ]$V7

Epic2_NarrowPeaks <- Epic2_peaks[width(Epic2_peaks) <= 2000, ]
Epic2_BroadPeaks <- Epic2_peaks[width(Epic2_peaks) > 2000, ]

#### Make a random control region file
cmd <- paste0("bedtools shuffle -i ", "DataFiles/Epic2_peaks_sorted.bed",
  " -g ", "DataFiles/mm39.chrom.sizes",
  " -excl ", "DataFiles/Epic2_peaks_sorted.bed",
  " -chrom -noOverlapping -seed 42 | sort -k1,1 -k2,2n > ", "DataFiles/NoEpic2PeakRegion.bed")
system(cmd)

cmd <- paste0("bedtools map -c 4 -o mean -a ",
              "DataFiles/NoEpic2PeakRegion.bed",
              " -b ", "DataFiles/H3K9me3_IgG_subtracted.bedgraph",
              " > ", "DataFiles/NoEpic2Region_methylSignal.bed")
system(cmd)
MethylSignal_noPeakRegion <- read.table("DataFiles/NoEpic2Region_methylSignal.bed") %>% rename(MeanIntensity = V7)
mean(MethylSignal_noPeakRegion$MeanIntensity)

#### Plot the mean methylation signal over these regions
narrow <- as.data.frame(Epic2_NarrowPeaks) %>% select(MeanIntensity) %>% mutate(type = "Narrow")
broad <- as.data.frame(Epic2_BroadPeaks) %>% select(MeanIntensity) %>% mutate(type = 'Broad')
NoPeak <- MethylSignal_noPeakRegion %>% select(MeanIntensity) %>% mutate(type = 'NoPeak')
df_to_plot <- rbind(narrow, broad, NoPeak)

ggplot(df_to_plot, aes(x = type, y = MeanIntensity)) +
  geom_bar(stat = 'summary', fun = 'mean', color = 'black', linewidth = 1, position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, size = 1, position = position_dodge(width = 0.8)) +
  labs(x = NULL)

ggplot(df_to_plot, aes(x = type, y = MeanIntensity)) +
  geom_boxplot(outliers = FALSE, linewidth = 1, fill = 'darkolivegreen3') +
  labs(x = NULL)

t.test(narrow$MeanIntensity, broad$MeanIntensity)

rm(narrow, broad, NoPeak)
###### PART II: Analyzing genes near different types of genes #######
### Categorized geens based on whether they overlap narrow peaks, broad domains or both
all_genes_10kb_extended <- resize(all_genes, width = width(all_genes) + 20000, fix = 'center')

all_genes_10kb_extended$NarrowPeakNum <- countOverlaps(all_genes_10kb_extended, Epic2_NarrowPeaks)
all_genes_10kb_extended$BroadPeakNum <- countOverlaps(all_genes_10kb_extended, Epic2_BroadPeaks)

### Type 1: genes only overlapping narrow peaks
genes_Nr <- subset(all_genes_10kb_extended, NarrowPeakNum > 0 & BroadPeakNum == 0)
### Type 2: genes only overlapping broad peaks
genes_Br <- subset(all_genes_10kb_extended, NarrowPeakNum == 0 & BroadPeakNum > 0)
### Type 3: genes only overlapping both peak types
genes_NrBr <- subset(all_genes_10kb_extended, NarrowPeakNum > 0 & BroadPeakNum > 0)

####### GO analysis on three gene types
library(clusterProfiler)
library(org.Mm.eg.db)

GO_list <- list()
GO_list[['Narrow']] <- enrichGO(genes_Nr$Symbol, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
GO_list[['Broad']] <- enrichGO(genes_Br$Symbol, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
GO_list[['Both']] <- enrichGO(genes_NrBr$Symbol, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)

write.xlsx(GO_list[['Narrow']], "DataFiles/GOAnalyses/GO_Genes_w_onlyNarrowPeaks.xlsx", rowNames = TRUE)
write.xlsx(GO_list[['Broad']], "DataFiles/GOAnalyses/GO_Genes_w_onlyBroadPeaks.xlsx", rowNames = TRUE)
write.xlsx(GO_list[['Both']], "DataFiles/GOAnalyses/GO_Genes_w_BothPeaks.xlsx", rowNames = TRUE)

barplot(GO_list[['Narrow']])
barplot(GO_list[['Broad']])
barplot(GO_list[['Both']])

####### Consider only narrow peaks regardless broad peaks
genes_Nr <- subset(all_genes_10kb_extended, NarrowPeakNum > 0)
length(genes_Nr)
go <- enrichGO(genes_Nr$Symbol, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go, showCategory = 10)
write.xlsx(as.data.frame(go), "DataFiles/GOAnalyses/GO_Genes_w_NarrowPeaks.xlsx", rowNames = TRUE)

#### Make barplot ####
## For genes with narrow peaks
selectedGO <- c('GO:0034329', 'GO:0050807', 'GO:0016358', 'GO:0097485', 'GO:0031589',
                'GO:0030010', 'GO:0043616', 'GO:0032253', 'GO:0043062', 'GO:0008544')

go_table <- read.xlsx("DataFiles/GOAnalyses/GO_Genes_w_NarrowPeaks.xlsx")
df_to_plot <- go_table %>% arrange(Count) %>% mutate(Description = factor(Description, levels = Description)) %>%
  mutate(LogPadjust = -log10(p.adjust)) %>%
  subset(ID %in% selectedGO)

p1 <- ggplot(df_to_plot, aes(x = Description, y = Count, fill = LogPadjust)) +
  geom_bar(stat = 'identity', color = 'black') +
  coord_flip() +
  scale_fill_gradient(low = 'white', high = 'darkgreen', limits = c(0, 9)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(x = NULL, y = "Gene count", fill = "-log10(p.adjust)", title = 'Genes w/ narrow peaks') +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12))

## For genes with only broad peaks
selectedGO <- c('GO:0019236', 'GO:0006805', 'GO:0097524', 'GO:0016229', 'GO:0016757')

go_table <- read.xlsx("DataFiles/GOAnalyses/GO_Genes_w_onlyBroadPeaks.xlsx")
df_to_plot <- go_table %>% arrange(Count) %>% mutate(Description = factor(Description, levels = Description)) %>%
  mutate(LogPadjust = -log10(p.adjust)) %>%
  subset(ID %in% selectedGO)

p2 <- ggplot(df_to_plot, aes(x = Description, y = Count, fill = LogPadjust)) +
  geom_bar(stat = 'identity', color = 'black') +
  coord_flip() +
  scale_fill_gradient(low = 'white', high = 'darkgreen', limits = c(0, 9), na.value = 'darkgreen') +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(x = NULL, y = "Gene count", fill = "-log10(p.adjust)", title = 'Genes w/ only broad peaks') +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = 'none')

p1/p2 + plot_layout(heights = c(2,1))

####### Check size distribution of peaks falling within 10kb of genes  #######
## Decide definition of narrow peaks based on this distribution
all_genes_10kb_extended <- resize(all_genes, width = width(all_genes) + 20000, fix = 'center')
gene_associated_peaks <- subset(Epic2_peaks, countOverlaps(Epic2_peaks, all_genes_10kb_extended) > 0)

summary(width(Epic2_peaks))

df_to_plot <- data.frame(Size = width(gene_associated_peaks), Score = gene_associated_peaks$Score)
df_to_plot$Bins <- cut(df_to_plot$Size, breaks = c(-Inf, 1000, 2000, 5000, 10000, 50000, Inf),
                       labels = c('< 1kb', '1-2kb', '2-5kb', '5-10kb', '10-50kb', '> 50kb'))
histogram <- ggplot(df_to_plot, aes(x = Bins)) +
  geom_bar(fill = 'black') +
  labs(x = 'Peak Length', y = 'Count')

table_df <- as.data.frame(t(fivenum(width(gene_associated_peaks))))
colnames(table_df) <- c('Min', 'Q1', 'Median', 'Q3', 'Max')
table_plot <- tableGrob(table_df, rows = NULL)
grid.arrange(table_plot, histogram, ncol = 1, heights = c(1, 4))

rm(table_df, table_plot, histogram)


####### Check RNA-seq fold change for three gene types  #######
RNAseq_E16_dataAll <- read.xlsx('/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/RNAseq_E16_genotypeOnly_noH614_DEgenes_annotated.xlsx', 
                                sheet = 1)

RNAseq_E16_dataAll <- RNAseq_E16_dataAll %>%
  mutate(PeakType = case_when(GeneID %in% genes_Nr$Symbol ~ 'Narrow',
                              GeneID %in% genes_Br$Symbol ~ 'Broad',
                              GeneID %in% genes_NrBr$Symbol ~ 'Both',
                              .default = 'No Peak')) %>%
  mutate(Tag = if_else(padj > 0.05 | abs(log2FoldChange) < 1, 'NS', if_else(log2FoldChange > 1, 'Up', 'Down'))) %>%
  mutate(Tag = factor(Tag, levels = c('NS', 'Down', 'Up')))

df_to_plot <- subset(RNAseq_E16_dataAll, padj > 1e-100)
ggplot(df_to_plot, aes(x = -log10(padj), y = log2FoldChange, colour = Tag)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c('gray50', 'blue', 'red')) +
  scale_x_continuous(limits = c(0, 40)) +
  facet_wrap(~PeakType, nrow = 1, axes = 'all')

df_to_plot_2 <- subset(RNAseq_E16_dataAll, GeneID %in% genes_Domain10kb$Symbol)
ggplot(df_to_plot_2, aes(x = -log10(padj), y = log2FoldChange, colour = Tag)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c('gray50', 'blue', 'red')) +
  scale_x_continuous(limits = c(0, 40))

