library(tidyverse)
library(openxlsx)

## For dReg peaks, take 300bp region from each dReg peak (150 upstream and 150 downstream from peak center)
#This step is already done in another script

## Prepare gene IDs for promoter region motif analysis
# Just a text file with a single column of ENSEMBL IDs
RNAseq_E16_DEs <- read.xlsx('/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/RNAseq_E16_genotypeOnly_noH614_DEgenes_annotated.xlsx', 
                            sheet = 2)
RNAseq_UpGenes <- subset(RNAseq_E16_DEs, log2FoldChange > 1, select = ensembl, drop = TRUE)

write.table(RNAseq_UpGenes, file = 'MotifAnalysisFiles/RNAseq_E16_upGenes.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)


PreciseSeqFC_table <- read.xlsx("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/Precise2024_FoldChangeTable_annotated.xlsx")
PreciseSeqFC_table$Category <- factor(PreciseSeqFC_table$Category, levels = c('upPI_nsGB', 'nsPI_upGB', 'downPI_upGB', 'downPI_nsGB', 'others'))
Nascent_GOI <- list()
for (name in unique(PreciseSeqFC_table$Category)){
  Nascent_GOI[[name]] <- subset(PreciseSeqFC_table, subset = Category == name, select = GeneID, drop = TRUE)
  write.table(Nascent_GOI[[name]], file = paste0('MotifAnalysisFiles/PreciseSeq_', name, '.txt'), row.names = FALSE, col.names = FALSE, quote = FALSE)
}


## Get the closest peak from different gene categories (especially upregulated genes)
#dReg peaks
dReg_distance <- read.table("../dREG/Closest_KOexclu_dRegPeaktoPromoters.bed", sep = "\t") %>%
  select(c(V4, V7:V15)) %>%
  rename(GeneID = V4, Symbol = V7, Biotype = V8, PeakScore = V13, Distance = V15,
         seqnames = V9, start = V10, end = V11, peakID = V12, strand = V14) %>%
  filter(PeakScore > 0)

#Epic2 peaks
Epic2_peaks <- read.table("/Users/kb744/Tumbar Lab Local/CUTandRUN/Functional_Analysis/DataFiles/Methyl_Annotation/ClosestEpic2Peak_toPromoters.bed") %>%
  select(c(V4, V7:V13)) %>%
  rename(GeneID = V4, PeakScore = V11, Distance = V13,
         seqnames = V7, start = V8, end = V9, peakID = V10, strand = V12) %>%
  filter(PeakScore > 0)


CB_export_cloest_peak <- function(peak, gList, fileName){
  dt <- subset(peak, GeneID %in% gList)
  dt <- dt[, c('seqnames', 'start', 'end', 'GeneID', 'PeakScore', 'strand')]
  write.table(dt, file = fileName, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
rm(CB_export_cloest_dReg)

CB_export_cloest_peak(dReg_distance, RNAseq_UpGenes, fileName = 'MotifAnalysisFiles/dRegPeak_KOexclu_RNAseqUp_genes.bed')
CB_export_cloest_peak(Epic2_peaks, RNAseq_UpGenes, fileName = 'MotifAnalysisFiles/H3K9me3Peak_RNAseqUp_genes.bed')


