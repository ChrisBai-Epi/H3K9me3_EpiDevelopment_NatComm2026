library(ChIPseeker)
library(tidyverse)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)

#Read in peak files
Epic2_peaks <- read.table("DataFiles/epic2_peaks_H3K9me3_E15.bed", header = TRUE)
Epic2_peaks$PeakID <- paste0("Epic", row.names(Epic2_peaks))
Epic2_peaks <- Epic2_peaks[, c(1, 2, 3, 11, 5, 6, 7, 8, 9, 10, 4)]
Epic2_peaks <- GRanges(Epic2_peaks)

#Meta of Epic2 peaks
library(grid)
library(gridExtra)

length(Epic2_peaks)
summary(width(Epic2_peaks))

df_to_plot <- data.frame(Size = width(Epic2_peaks), Score = Epic2_peaks$Score)
df_to_plot$Bins <- cut(df_to_plot$Size, breaks = c(-Inf, 1000, 2000, 5000, 10000, 50000, Inf),
                       labels = c('< 1kb', '1-2kb', '2-5kb', '5-10kb', '10-50kb', '> 50kb'))
histogram <- ggplot(df_to_plot, aes(x = Bins)) +
  geom_bar(fill = 'black') +
  labs(x = 'Peak Length', y = 'Count')

table_df <- as.data.frame(t(fivenum(width(Epic2_peaks))))
colnames(table_df) <- c('Min', 'Q1', 'Median', 'Q3', 'Max')
table_plot <- tableGrob(table_df, rows = NULL, base_family = 'Arial')
grid.arrange(table_plot, histogram, ncol = 1, heights = c(1, 4))

rm(table_df, table_plot, histogram)


#Peak annotation
peakAnno <- annotatePeak(Epic2_peaks, tssRegion = c(-500, 500), TxDb = TxDb.Mmusculus.UCSC.mm39.knownGene, annoDb = "org.Mm.eg.db",
                         genomicAnnotationPriority = c("Promoter", "Exon", "Intron", "5UTR", "3UTR", "Downstream", "Intergenic"))
#plotAnnoPie(peakAnno, cex = 0.7)

AnnoSummary <- data.frame(Feature = c("Promoter", "Exon", "Intron", "Intergenic"),
                          Percentage_raw = c(2.3241, 6.5154, 23.3578, 67.8027))
AnnoSummary$Percentage <- round(AnnoSummary$Percentage_raw, digits = 1)
ggplot(AnnoSummary, aes(x='', y = Percentage, fill = Feature)) +
  geom_bar(stat = 'identity', width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Accent") +
  geom_text(aes(label = paste0(Feature, '\n', Percentage)), 
            position = position_stack(vjust = 0.5), 
            color = "black", size = 4)+
  theme_void() +
  theme(legend.position = 'none')

#plotCorrelation of BAM files
system("deeptools --version")

cmd <- paste0('plotCorrelation -in ', '../QualityControl/multiBamsummary_all.npz',
              ' -c pearson -p heatmap --plotNumbers --removeOutliers -o ',
              '../QualityControl/Heatmap_Correlation_CUTandRUN_2.pdf',
              ' --colorMap RdBu_r')
system(cmd)

cmd <- paste0('plotCorrelation -in ', '../QualityControl/multiBamsummary_all.npz',
              ' -c spearman -p heatmap --plotNumbers --removeOutliers -o ',
              '../QualityControl/Heatmap_Correlation_CUTandRUN.pdf',
              ' --colorMap RdBu_r')
system(cmd)


