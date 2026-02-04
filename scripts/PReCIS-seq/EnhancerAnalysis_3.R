library(tidyverse)
library(openxlsx)
library(ggrastr)
library(GenomicRanges)

## This script focuses on computing for the methylation signal and nascent transcription signal
## around dREGs

## Score cutoff set at 0.4
nonPr_dREG_peaks <- lapply(dRegPeaks, function(gr){
  pr <- promoters(AllGenes, upstream = 500, downstream = 500)
  overlapped <- findOverlaps(gr, pr)
  np_gr <- gr[setdiff(seq_along(gr), queryHits(overlapped))]
  np_gr_high <- subset(np_gr, score > 0.4)
})

nonPr_dREG_peaks$CT_IP$exclusive <- TRUE
nonPr_dREG_peaks$KO_IP$exclusive <- TRUE
idx <- findOverlaps(nonPr_dREG_peaks$CT_IP, nonPr_dREG_peaks$KO_IP)
nonPr_dREG_peaks$CT_IP$exclusive[unique(queryHits(idx))] <- FALSE
nonPr_dREG_peaks$KO_IP$exclusive[unique(subjectHits(idx))] <- FALSE

nonPr_dREG_peaks$CT_IP$tag <- 'Ctrl'
nonPr_dREG_peaks$KO_IP$tag <- 'TKO'

CB_export_dREG(nonPr_dREG_peaks$CT_IP, "../dREG/NonPr_dREGs_CT_IP_1kb.bed", width = 1000)
CB_export_dREG(nonPr_dREG_peaks$KO_IP, "../dREG/NonPr_dREGs_KO_IP_1kb.bed", width = 1000)

CB_export_dREG(nonPr_dREG_peaks$CT_IP, "../dREG/NonPr_dREGs_CT_IP_2kb.bed", width = 2000)
CB_export_dREG(nonPr_dREG_peaks$KO_IP, "../dREG/NonPr_dREGs_KO_IP_2kb.bed", width = 2000)


## Using bedtools map to calculate H3K9me3 signal signal on dREG regions
cmd <- paste0("bedtools map -c 4 -o sum -a ",
              "../dREG/NonPr_dREGs_CT_IP_2kb.bed",
              " -b ", "../BrowserTracks/H3K9me3_IgG_subtracted.bedgraph",
              " > ", "../dREG/NonPr_dREGs_CT_IP_2kb_methylinfo.bed")
system(cmd)

cmd <- paste0("bedtools map -c 4 -o sum -a ",
              "../dREG/NonPr_dREGs_KO_IP_2kb.bed",
              " -b ", "../BrowserTracks/H3K9me3_IgG_subtracted.bedgraph",
              " > ", "../dREG/NonPr_dREGs_KO_IP_2kb_methylinfo.bed")
system(cmd)

## Add methyl info to nonPr_dREG_peak objects
methylInfo <- read.table("../dREG/NonPr_dREGs_CT_IP_2kb_methylinfo.bed") %>%
  rename(methyl_signal = V9) %>% column_to_rownames(var = 'V4')
nonPr_dREG_peaks$CT_IP$methyl_signal <- methylInfo[nonPr_dREG_peaks$CT_IP$name, 'methyl_signal']

methylInfo <- read.table("../dREG/NonPr_dREGs_KO_IP_2kb_methylinfo.bed") %>%
  rename(methyl_signal = V9) %>% column_to_rownames(var = 'V4')
nonPr_dREG_peaks$KO_IP$methyl_signal <- methylInfo[nonPr_dREG_peaks$KO_IP$name, 'methyl_signal']

set.seed(428)
random_CT_dREGs <- nonPr_dREG_peaks$CT_IP[sample(length(nonPr_dREG_peaks$CT_IP), 2500)]
random_CT_dREGs$tag <- 'Random_CT'
KO_exclusive_dRegPeaks <- subset(nonPr_dREG_peaks$KO_IP, exclusive)
KO_exclusive_dRegPeaks$tag <- 'KO_exclusive'
top_KO_dREGs <- head(nonPr_dREG_peaks$KO_IP[order(nonPr_dREG_peaks$KO_IP$methyl_signal, decreasing = TRUE)], 500)
top_KO_dREGs <- head(KO_exclusive_dRegPeaks[order(KO_exclusive_dRegPeaks$methyl_signal, decreasing = TRUE)], 500)
top_KO_dREGs$tag <- 'Top_methylated'

df_to_plot <- as.data.frame(c(random_CT_dREGs, KO_exclusive_dRegPeaks, top_KO_dREGs)) %>%
  mutate(tag = factor(tag, levels = c('Random_CT', 'KO_exclusive', 'Top_methylated')))
df_to_plot <- subset(df_to_plot, tag != 'KO_exclusive') %>%
  mutate(tag = recode(tag, "Random_CT" = "Ctrl dREGs", "Top_methylated" = "TKO only dREGs"))
ggplot(df_to_plot, aes(x = tag, y = methyl_signal, fill = tag)) +
  geom_bar(stat = 'summary', fun = 'mean', width = 0.8, color = 'black', size = 1) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, size =1) +
  scale_fill_manual(values = genoColor) +
  geom_hline(yintercept = 0, color = 'black', linewidth = 1, linetype = 'dashed') +
  labs(x = NULL, y = 'Methylation Signal') +
  theme(legend.position = 'none')

rm(random_CT_dREGs, KO_exclusive_dRegPeaks, top_KO_dREGs, methylInfo)

CT_dREGs_copy <- nonPr_dREG_peaks$CT_IP[sample(length(nonPr_dREG_peaks$CT_IP), 2500)]
CT_dREGs_copy$tag <- 'Ctrl'
CT_dREGs_copy <-head(CT_dREGs_copy[order(CT_dREGs_copy$methyl_signal, decreasing = TRUE)], 2500)
CT_dREGs_copy$rank <- seq_along(CT_dREGs_copy)

KO_exclusive_dREGs <- subset(nonPr_dREG_peaks$KO_IP, exclusive)
KO_exclusive_dREGs <- head(KO_exclusive_dREGs[order(KO_exclusive_dREGs$methyl_signal, decreasing = TRUE)], 2500)
KO_exclusive_dREGs$tag <- 'TKO'
KO_exclusive_dREGs$rank <- seq_along(KO_exclusive_dREGs)

df_to_plot <- as.data.frame(c(CT_dREGs_copy, KO_exclusive_dREGs))
df_to_plot$Bins <- cut(df_to_plot$rank, breaks = c(0, 500, 1000, 1500, 2000, Inf),
                       labels = c('1-500', '500-1000', '1000-1500', '1500-2000', '>2000'))

ggplot(df_to_plot, aes(x = Bins, y = methyl_signal, fill = tag)) +
  geom_bar(stat = 'summary', fun = 'mean', width = 0.8, color = 'black', size = 1, position = position_dodge(width = 0.8)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, size = 1, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = genoColor) +
  labs(x = 'Rank', y = 'Methylation signal', fill = 'Group')

rm(CT_dREGs_copy, KO_exclusive_dREGs)

## Barplot of methylated dREGs vs unmethylated dREGs; obtaining process is in Enhancer_plotting_dREGs.R
methylated_dREGs$tag <- 'Methylated'
Unmethylated_dREGs$tag <- 'Unmethylated'
CT_dREGs_copy$tag <- 'Ctrl dREGs'
df_to_plot <- as.data.frame(c(methylated_dREGs, Unmethylated_dREGs, CT_dREGs_copy))
df_to_plot$tag <- factor(df_to_plot$tag, levels = c('Ctrl dREGs', 'Unmethylated', 'Methylated'))

ggplot(df_to_plot, aes(x = tag, y = methyl_signal, fill = tag)) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = c('darkgreen', 'lightyellow', 'gold')) +
  labs(x = NULL, y = 'H3K9me3 signal', fill = 'Group')

t.test(methylated_dREGs$methyl_signal, CT_dREGs_copy$methyl_signal)
t.test(Unmethylated_dREGs$methyl_signal, CT_dREGs_copy$methyl_signal)

ggplot(df_to_plot, aes(x = tag, y = methyl_signal, fill = tag)) +
  geom_bar(stat = 'summary', fun = 'mean', width = 0.8, color = 'black', size = 1, position = position_dodge(width = 0.8)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, size = 1, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c('darkgreen', 'lightyellow', 'gold')) +
  labs(x = NULL, y = 'H3K9me3 signal', fill = 'Group')

## Using bedtools map to calculate nascent transcription signal on dREG regions
cmd <- paste0("bedtools map -c 4 -o sum -a ",
              "../dREG/NonPr_dREGs_CT_IP_1kb.bed",
              " -b ", "../BrowserTracks/Bedgraphs/CT_IP_fwd.bedgraph",
              " > ", "../dREG/NonPr_dREGs_CT_IP_1kb_nascent_fwd.bed")
system(cmd)
cmd <- paste0("bedtools map -c 4 -o sum -a ",
              "../dREG/NonPr_dREGs_CT_IP_1kb.bed",
              " -b ", "../BrowserTracks/Bedgraphs/CT_IP_rev.bedgraph",
              " > ", "../dREG/NonPr_dREGs_CT_IP_1kb_nascent_rev.bed")
system(cmd)

cmd <- paste0("bedtools map -c 4 -o sum -a ",
              "../dREG/NonPr_dREGs_KO_IP_1kb.bed",
              " -b ", "../BrowserTracks/Bedgraphs/KO_IP_fwd.bedgraph",
              " > ", "../dREG/NonPr_dREGs_KO_IP_1kb_nascent_fwd.bed")
system(cmd)
cmd <- paste0("bedtools map -c 4 -o sum -a ",
              "../dREG/NonPr_dREGs_KO_IP_1kb.bed",
              " -b ", "../BrowserTracks/Bedgraphs/KO_IP_rev.bedgraph",
              " > ", "../dREG/NonPr_dREGs_KO_IP_1kb_nascent_rev.bed")
system(cmd)

## Add Nascent transcription info to nonPr_dREG_peaks object
nascent_info <- read.table("../dREG/NonPr_dREGs_CT_IP_1kb_nascent_fwd.bed") %>%
  rename(signal = V9) %>% column_to_rownames(var = 'V4')
nonPr_dREG_peaks$CT_IP$nascent_fwd <- nascent_info[nonPr_dREG_peaks$CT_IP$name, 'signal']

nascent_info <- read.table("../dREG/NonPr_dREGs_CT_IP_1kb_nascent_rev.bed") %>%
  rename(signal = V9) %>% column_to_rownames(var = 'V4')
nonPr_dREG_peaks$CT_IP$nascent_rev <- nascent_info[nonPr_dREG_peaks$CT_IP$name, 'signal']

nascent_info <- read.table("../dREG/NonPr_dREGs_KO_IP_1kb_nascent_fwd.bed") %>%
  rename(signal = V9) %>% column_to_rownames(var = 'V4')
nonPr_dREG_peaks$KO_IP$nascent_fwd <- nascent_info[nonPr_dREG_peaks$KO_IP$name, 'signal']

nascent_info <- read.table("../dREG/NonPr_dREGs_KO_IP_1kb_nascent_rev.bed") %>%
  rename(signal = V9) %>% column_to_rownames(var = 'V4')
nonPr_dREG_peaks$KO_IP$nascent_rev <- nascent_info[nonPr_dREG_peaks$KO_IP$name, 'signal']

nonPr_dREG_peaks <- lapply(nonPr_dREG_peaks, function(gr){
  gr$dREG_signal <- gr$nascent_fwd - gr$nascent_rev
  return(gr)
})

top_CT_dREGs <- head(nonPr_dREG_peaks$CT_IP[order(nonPr_dREG_peaks$CT_IP$score, decreasing = TRUE)], 1000)
CB_export_dREG(top_CT_dREGs, "../dREG/NonPr_dREGs_CT_IP_topScore.bed")

## Plot nascent transcription for dREGs
command <- paste0("computeMatrix reference-point -S ", 
                  "../BrowserTracks/SpikeIn_Norm_BigWigs/CT_IP_fwd.bw", " ",
                  "../BrowserTracks/SpikeIn_Norm_BigWigs/CT_IP_rev.bw",
                  " -R ", "../dREG/NonPr_dREGs_CT_IP_topScore.bed",
                  " -o ", "../dREG/Heatmap/NascentRNA_matrix_dREGs.gz",
                  " --referencePoint center --binSize 10 --sortRegions descend -p 15 -a 500 -b 500")
system(command)

command <- paste0("plotHeatmap -m ", "../dREG/Heatmap/NascentRNA_matrix_dREGs.gz",
                  " -o ", "../dREG/Heatmap/NascentRNA_dREGs.pdf",
                  " --colorMap RdBu_r --samplesLabel H3K9me3 --samplesLabel Forward Reverse --regionsLabel dREGs --refPointLabel PeakCenter",
                  " --heatmapWidth 8 --heatmapHeight 16 --averageTypeSummaryPlot mean --perGroup --zMax 1.5 --zMin -1.5")
system(command)

## Plot nascent transcription for KO exclusive dREGs
KO_exclusive_dREGs <- subset(nonPr_dREG_peaks$KO_IP, exclusive)
CB_export_dREG(KO_exclusive_dREGs, "../dREG/NonPr_dREGs_KO_IP_exclusivePeaks.bed")

command <- paste0("computeMatrix reference-point -S ", 
                  "../BrowserTracks/SpikeIn_Norm_BigWigs/KO_IP_fwd.bw", " ",
                  "../BrowserTracks/SpikeIn_Norm_BigWigs/KO_IP_rev.bw",
                  " -R ", "../dREG/NonPr_dREGs_KO_IP_exclusivePeaks.bed",
                  " -o ", "../dREG/Heatmap/NascentRNA_matrix_KO_dREGs.gz",
                  " --referencePoint center --binSize 10 --sortRegions descend -p 15 -a 500 -b 500")
system(command)

command <- paste0("plotHeatmap -m ", "../dREG/Heatmap/NascentRNA_matrix_KO_dREGs.gz",
                  " -o ", "../dREG/Heatmap/NascentRNA_KOex_dREGs.pdf",
                  " --colorMap RdBu_r --samplesLabel H3K9me3 --samplesLabel Forward Reverse --regionsLabel dREGs --refPointLabel PeakCenter",
                  " --heatmapWidth 8 --heatmapHeight 16 --averageTypeSummaryPlot mean --perGroup --zMax 1.5 --zMin -1.5")
system(command)



