library(tidyverse)
library(GenomicRanges)

#### Obtain Epic2_peaks (separate as GRange object) with QCandProcessing.R
Epic2_peaks[['P12']] <- NULL

#### Meta of Epic2 peaks, get numbers here, and create tables with InDesign
Epic2_peaks_all <- c(Epic2_peaks[['E12']], Epic2_peaks[['E14']], Epic2_peaks[['E16']])
length(Epic2_peaks_all)
summary(width(Epic2_peaks[['E12']]))
summary(width(Epic2_peaks[['E14']]))
summary(width(Epic2_peaks[['E16']]))

#### Make histograms
df_to_plot <- data.frame(ID = Epic2_peaks_all$peakID,
                         Size = width(Epic2_peaks_all),
                         Score = Epic2_peaks_all$score)
df_to_plot <- df_to_plot %>%
  mutate(Bins = cut(Size, breaks = c(-Inf, 2000, 5000, 10000, 50000, Inf),
                    labels = c('< 2kb', '2-5kb', '5-10kb', '10-50kb', '> 50kb'))) %>%
  mutate(Stage = str_extract(ID, '^E\\d+')) %>%
  group_by(Stage, Bins) %>% summarise(Count = n()) %>% mutate(Frac = Count/sum(Count)) %>% ungroup()

ggplot(df_to_plot, aes(x = Bins, y = Frac)) +
  geom_bar(fill = 'black', stat = 'identity') +
  facet_wrap(~Stage) +
  labs(x = 'Peak Length', y = 'Count')

#### Genomic Annotation of peaks
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(org.Mm.eg.db)

peakAnno <- annotatePeak(Epic2_peaks_all, tssRegion = c(-500, 500), TxDb = TxDb.Mmusculus.UCSC.mm39.knownGene, annoDb = "org.Mm.eg.db",
                         genomicAnnotationPriority = c("Promoter", "Exon", "Intron", "5UTR", "3UTR", "Downstream", "Intergenic"))
peakAnno@annoStat

#### The percentage number here is hard-coded with the values from peakAnno@annoStat
AnnoSummary <- data.frame(Feature = c("Promoter", "Exon", "Intron", "Intergenic"),
                          Percentage_raw = c(8.75, 17.43, 25.88, 47.94))
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

rm(peakAnno, AnnoSummary)

#### Plot types of cCREs overlapping Epic2 peaks (methylated)
cCREs <- read.table("DataFiles/mm39_cCREs_all_lifeovered.bed")
colnames(cCREs) <- c('seqnames', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'label')
cCREs <- GRanges(cCREs)
cCREs$label <- sub(",.*", "", cCREs$label)

cCREs$Methylated <- countOverlaps(cCREs, Epic2_peaks_all) > 0
Epic2_peaks_all$cCRE <- countOverlaps(Epic2_peaks_all, cCREs)
sum(Epic2_peaks_all$cCRE > 0)

df_to_plot <- as.data.frame(cCREs) %>% subset(Methylated) %>%
  group_by(label) %>% summarise(Count = n()) %>%
  mutate(Pct = Count / sum(Count))

ggplot(df_to_plot, aes(x = '', y = Pct, fill = label)) +
  geom_bar(stat = 'identity') +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_fill_brewer(palette = 'Set3') +
  labs(x = NULL, y = "Percentage")

rm(Epic2_peaks_all)

#### Check if dynamic peaks are more narrow ####
df_to_plot <- Epic2_union_peak %>% mutate(Width = abs(start - end))

ggplot(df_to_plot, aes(x = Dynamics, y = Width)) +
  geom_boxplot(outliers = FALSE)

