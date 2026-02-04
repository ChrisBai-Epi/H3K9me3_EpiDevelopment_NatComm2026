library(tidyverse)

## Obtain the top 500 dREGs for plotting
CT_dREGs_copy <- nonPr_dREG_peaks$CT_IP[sample(length(nonPr_dREG_peaks$CT_IP), 2500)]
CT_dREGs_copy$tag <- 'Ctrl'
CT_dREGs_copy <-head(CT_dREGs_copy[order(CT_dREGs_copy$methyl_signal, decreasing = TRUE)], 2500)
CT_dREGs_copy$rank <- seq_along(CT_dREGs_copy)

KO_exclusive_dREGs <- subset(nonPr_dREG_peaks$KO_IP, exclusive)
KO_exclusive_dREGs <- head(KO_exclusive_dREGs[order(KO_exclusive_dREGs$methyl_signal, decreasing = TRUE)], 2500)
KO_exclusive_dREGs$tag <- 'TKO'
KO_exclusive_dREGs$rank <- seq_along(KO_exclusive_dREGs)

CB_export_dREG(head(CT_dREGs_copy, 500), "../dREG/Heatmap_Apr2025/CT_dREGs_500.bed")
CB_export_dREG(head(KO_exclusive_dREGs, 500), "../dREG/Heatmap_Apr2025/KO_dREGs_500.bed")


## Plotting methylation signals
system2("deeptools", args = "--version")
command <- paste0("computeMatrix reference-point -S ", "../dREG/Heatmap/H3K9me3_IgG_subtracted.bw",
                  " -R ", "../dREG/Heatmap_Apr2025/CT_dREGs_500.bed", " ", "../dREG/Heatmap_Apr2025/KO_dREGs_500.bed",
                  " -o ", "../dREG/Heatmap_Apr2025/matrix_dREGs.gz",
                  " --referencePoint center --binSize 10 --sortRegions descend -p 15 -a 500 -b 500",
                  " --maxThreshold 10")
system(command)

system("gunzip ../dREG/Heatmap_Apr2025/matrix_dREGs.gz")
system("gzip ../dREG/Heatmap_Apr2025/matrix_dREGs")

command <- paste0("plotHeatmap -m ", "../dREG/Heatmap_Apr2025/matrix_dREGs.gz",
                  " -o ", "../dREG/Heatmap_Apr2025/Methylation_dREGs.pdf",
                  " --colorMap RdYlBu_r --samplesLabel H3K9me3 --regionsLabel Ctrl TKO --refPointLabel PeakCenter",
                  " --heatmapWidth 8 --heatmapHeight 16 --averageTypeSummaryPlot mean")
system(command)

##### Select methylation-regulated dREGs for Homer motif analysis
## Testing each individual KO values against Ctrl distribution
Ctrl_signals <- CT_dREGs_copy$methyl_signal[-1]
z_scores <- (KO_exclusive_dREGs$methyl_signal - mean(Ctrl_signals)) / sd(Ctrl_signals)
pvalues <- 1 - pnorm(z_scores)
significant <- pvalues < 0.05
KO_exclusive_dREGs$MethylTag[significant] <- 'methylated'
KO_exclusive_dREGs$MethylTag[!significant] <- 'unmethylated'
methylated_dREGs <- KO_exclusive_dREGs[significant]
Unmethylated_dREGs <- KO_exclusive_dREGs[!significant]

rm(Ctrl_signals, z_scores, pvalues, significant)

CB_PrepareBed_forHomer(methylated_dREGs, "MotifAnalysisFiles/methylated_dREGs.bed")
CB_PrepareBed_forHomer(Unmethylated_dREGs, "MotifAnalysisFiles/Unmethylated_dREGs.bed")

rm(methylated_dREGs, Unmethylated_dREGs)

a <- as.data.frame(KO_exclusive_dREGs) %>% group_by(epic2, MethylTag) %>% summarise(Count = n()) %>% ungroup

### Combine plus and minus strand and plot transcription on dRegs
system('Deeptools --version')

# Create combined bigwig files
command <- paste0('bigwigCompare -b1 ',
                  'SpikeIn_Norm_BigWigs/CT_IP_fwd.bw', ' -b2 ', 'SpikeIn_Norm_BigWigs/CT_IP_rev.bw',
                  ' --operation add -bs 5 -p 20 -of bigwig -o ', 'SpikeIn_Norm_BigWigs/CT_IP_combined.bw')
system(command)
command <- paste0('bigwigCompare -b1 ',
                  'SpikeIn_Norm_BigWigs/KO_IP_fwd.bw', ' -b2 ', 'SpikeIn_Norm_BigWigs/KO_IP_rev.bw',
                  ' --operation add -bs 5 -p 20 -of bigwig -o ', 'SpikeIn_Norm_BigWigs/KO_IP_combined.bw')
system(command)

# Plotting
command <- paste0("computeMatrix reference-point -S ", 
                  "SpikeIn_Norm_BigWigs/CT_IP_combined.bw", " ", "SpikeIn_Norm_BigWigs/KO_IP_combined.bw",
                  " -R ", "../dREG/NonPr_dREGs_KO_IP_exclusivePeaks.bed",
                  " -o ", "../dREG/Heatmap/NascentRNA_matrix_dREGs.gz",
                  " --referencePoint center --binSize 20 --sortRegions descend -p 15 -a 500 -b 500")
system(command)

command <- paste0("plotHeatmap -m ", "../dREG/Heatmap/NascentRNA_matrix_dREGs.gz",
                  " -o ", "../dREG/Heatmap/NascentRNA_dREGs_strandCombined.pdf",
                  " --colorMap RdBu_r --samplesLabel Ctrl TKO --regionsLabel dREGs --refPointLabel PeakCenter",
                  " --heatmapWidth 8 --heatmapHeight 16 --averageTypeSummaryPlot mean --perGroup --zMax 1.5 --zMin -1.5")
system(command)
