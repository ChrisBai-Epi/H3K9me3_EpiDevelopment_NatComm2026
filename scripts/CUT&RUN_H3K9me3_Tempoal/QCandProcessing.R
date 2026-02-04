## plotCorrelation of BAM files
system("deeptools --version")

cmd <- paste0('plotCorrelation -in ', 'DataFiles/coverageMatrix_peakRegion.npz',
              ' -c pearson -p heatmap --plotNumbers -o ',
              'DataFils/Heatmap_Correlation_peakRegion.pdf',
              ' --colorMap RdBu_r --zMin 0.8')
system(cmd)

cmd <- paste0('plotCorrelation -in ', 'DataFiles/coverageMatrix_bins.npz',
              ' -c spearman -p heatmap --plotNumbers -o ',
              'DataFiles/Heatmap_Correlation_bins.pdf',
              ' --colorMap RdBu_r')
system(cmd)

cmd <- paste0('plotCorrelation -in ', 'DataFiles/CorrelationMatrix_AllStages.npz',
              ' -c spearman -p heatmap --plotNumbers -o ',
              'DataFiles/Heatmap_Correlation_allStages_bins.pdf',
              ' --colorMap RdBu_r')
system(cmd)

## Plot PCA
cmd <- paste0('plotPCA -in ', 'DataFiles/coverageMatrix_peakRegion.npz',
              ' --rowCenter --plotFileFormat pdf -o ',
              'DataFiles/PeakRegions_PCA.pdf')
system(cmd)


#Process Epic2 peak calling results
library(tidyverse)
library(GenomicRanges)
library(openxlsx)

#Read in peak files
file_names <- list("E12" = '../Epic2Peaks/epic2_peaks_H3K9me3_E12.bed',
                   "E14" = '../Epic2Peaks/epic2_peaks_H3K9me3_E14.bed',
                   "E16" = '../Epic2Peaks/epic2_peaks_H3K9me3_E16.bed',
                   "P12" = '../Epic2Peaks/epic2_peaks_H3K9me3_P12.bed')

Epic2_peaks <- lapply(file_names, function(fn){
  peak <- read.table(fn, header = FALSE)
  names(peak) <- c('chrom', 'start', 'end', 'pvalue', 'score', 'strand', 'H3K9_Count', 'IgG_Count', 'FDR', 'log2FC')
  peak$peakID <- paste0('peak_', row.names(peak))
  peak <- peak[, c(1, 2, 3, 11, 5, 6, 7, 8, 9, 10, 4)]
})

for (stage in c("E12", "E14", "E16", "P12")){
  Epic2_peaks[[stage]]$peakID <- paste0(stage, '_', Epic2_peaks[[stage]]$peakID)
}

for (stage in c("E12", "E14", "E16", "P12")){
  write.table(Epic2_peaks[[stage]][,1:6], file = paste0("../Epic2Peaks/ForIGV_Epic2_peaks_", stage, ".bed"),
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
}

Epic2_peaks <- lapply(Epic2_peaks, function(df){
  gr <- GRanges(df)
})

rm(Epic2_peaks, file_names, stage)

local({
  peak <- read.table('../Epic2Peaks/epic2_peaks_H3K9me3_E12_test.bed', header = FALSE)
  names(peak) <- c('chrom', 'start', 'end', 'pvalue', 'score', 'strand', 'H3K9_Count', 'IgG_Count', 'FDR', 'log2FC')
  peak$peakID <- paste0('p', row.names(peak))
  peak <- peak[, c(1, 2, 3, 11, 5, 6)]
  write.table(peak, file = paste0("../Epic2Peaks//ForIGV_Epic2_peaks_E12_test.bed"),
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
})


### Process combined peak regions and H3K9me3 counts
Epic2_union_peak <- read.table('../Epic2Peaks/Union_epic2_H3K9me3_peaks.bed', header = FALSE)
colnames(Epic2_union_peak) <- c('chrom', 'start', 'end', 'Ex_ID', 'score', 'pvalue')
rawCounts <- read.table('DataFiles/rawCounts_sorted.txt')

Epic2_union_peak$PeakID <- paste0('Epic2_Peak_', row.names(Epic2_union_peak))
Epic2_union_peak <- Epic2_union_peak[, c(1, 2, 3, 7, 5, 6, 4)]
write.table(Epic2_union_peak, file = paste0("../Epic2Peaks/ForIGV_Epic2_union_peaks.bed"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')

Counts <- cbind(Epic2_union_peak, rawCounts[,4:14])
colnames(Counts)[8:18] <- c('E12_1_IgG', 'E12_2_H3K9', 'E12_3_H3K9',
                                      'E14_1_IgG', 'E14_2_H3K9', 'E14_3_H3K9', 'E14_4_H3K9',
                                      'E16_1_H3K9', 'E16_2_H3K9', 'E16_3_IgG', 'E16_4_H3K9') 
rm(rawCounts)



