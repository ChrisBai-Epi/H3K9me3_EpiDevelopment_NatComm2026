library(tidyverse)
library(openxlsx)
library(ggrastr)
library(GenomicRanges)

## Load dREG related data with EnhancerAnalysis.R
## This script focusing on analyzing how dREGs overlaps H3K9me3 peak and cCREs

#Prepare cCREs first by extracting the first 6 columns for UCSC liftover tool
cCREs <- read.table("../BrowserTracks/mm10_cCREs_all_ENCFF904ZZH.bed") %>%
  select(c(V1:V6, V10))
write.table(cCREs, file = "../BrowserTracks/mm10_cCREs_all_V1toV6.bed", col.names = FALSE, row.names = FALSE, quote = FALSE)

cCREs <- read.table("../BrowserTracks/mm39_cCREs_all_lifeovered.bed")
colnames(cCREs) <- c('seqnames', 'start', 'end', 'name', 'score', 'strand', 'label')
cCREs$score <- 0
cCRE_color_info <- read.table("../BrowserTracks/mm10_cCREs_all_ENCFF904ZZH.bed") %>% select(c(V4, V7:V9)) %>%
  rename(name = V4, thickStart = V7, thickEnd = V8, itemRgb = V9)
cCREs <- left_join(cCREs, cCRE_color_info, by = 'name')
cCREs <- cCREs[, c(1:6, 8, 9, 10, 7)]
cCREs[, 7:8] <- cCREs[, 2:3] 
write.table(cCREs, file = "../BrowserTracks/mm39_cCREs_all_lifeovered.bed", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

cCREs <- read.table("../BrowserTracks/mm39_cCREs_all_lifeovered.bed")
colnames(cCREs) <- c('seqnames', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'label')
cCREs_gr <- GRanges(cCREs)
cCREs_tb <- cCREs
rm(cCREs)

Epic2_peaks <- read.table("../BrowserTracks/Epic2_peaks.bed")
colnames(Epic2_peaks) <- c('seqnames', 'start', 'end', 'name', 'score', 'strand')
Epic2_peaks <- GRanges(Epic2_peaks)

## Analysize the percentage of dREGs covered by H3K9me3
nonPr_dREG_peaks <- lapply(nonPr_dREG_peaks, function(gr){
  idx <- findOverlaps(gr, Epic2_peaks)
  gr$epic2 <- FALSE
  gr$epic2[unique(queryHits(idx))] <- TRUE
  return(gr)
})

df_to_plot <- as.data.frame(c(nonPr_dREG_peaks$CT_IP, nonPr_dREG_peaks$KO_IP))
df_to_plot <- df_to_plot %>% group_by(tag, exclusive, epic2) %>% summarise(Count = n()) %>%
  mutate(Frac = Count / sum(Count)) %>% subset(epic2) %>% mutate(Type = if_else(exclusive, 'Exclusive', 'Common'))

ggplot(df_to_plot, aes(x = Type, y = Frac, fill = tag)) +
  geom_bar(stat = 'identity', position = position_dodge(), color = 'black', size = 0.75) +
  scale_y_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
  scale_fill_manual(values = genoColor) +
  labs(x = NULL, y = 'Percent of methylated dREGs', fill = NULL)

## Analyzing how dREGs overlap cCREs
cCREs_tb$label <- sub(",.*", "", cCREs_tb$label)
cCREs_gr$label <- sub(",.*", "", cCREs_gr$label)

nonPr_dREG_peaks <- lapply(nonPr_dREG_peaks, function(gr){
  gr$cCRE <- findOverlaps(gr, cCREs_gr, select = 'first')
  gr$cCRE <- cCREs_gr$label[gr$cCRE]
  return(gr)
})

df_to_plot <- as.data.frame(nonPr_dREG_peaks$KO_IP)
df_to_plot <- df_to_plot %>% mutate(label = if_else(cCRE != 'pELS' & cCRE != 'dELS', 'others', cCRE)) %>%
  mutate(label = if_else(is.na(label), 'non-cCRE', label)) %>%
  mutate(tag = if_else(exclusive, 'Exclusive', 'Common')) %>%
  group_by(tag, label) %>% summarise(Count = n()) %>% mutate(Frac = Count / sum(Count)) %>% ungroup() %>%
  mutate(label = factor(label, levels = c('non-cCRE', 'others', 'pELS', 'dELS')))

ggplot(df_to_plot, aes(x = tag, y = Count, fill = label)) +
  geom_bar(stat = 'identity', color = 'black', size = 0.75) +
  scale_fill_brewer(palette = 'Set3') +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = 'Counts', fill = 'cCRE')

## Make circular barplot
df_to_plot <- as.data.frame(nonPr_dREG_peaks$KO_IP)
df_to_plot <- df_to_plot %>% mutate(label = if_else(cCRE != 'pELS' & cCRE != 'dELS', 'others', cCRE)) %>%
  mutate(label = if_else(is.na(label), 'non-cCRE', label)) %>%
  mutate(tag = if_else(exclusive, 'Exclusive', 'Common')) %>%
  group_by(tag, label) %>% summarise(Count = n()) %>% ungroup() %>%
  mutate(Frac = Count / sum(Count), ring = 1.5)

df_to_plot_inner <- df_to_plot %>% group_by(tag) %>% summarise(Count = sum(Count)) %>%
  mutate(Frac = Count / sum(Count)) %>% mutate (label = 'Total', ring = 0.5)
df_to_plot <- bind_rows(df_to_plot, df_to_plot_inner) %>%
  mutate(label_2 = factor(interaction(label, tag),
                          levels = c('Total.Common', 'Total.Exclusive',
                                     'dELS.Common', 'pELS.Common', 'others.Common', 'non-cCRE.Common',
                                     'dELS.Exclusive', 'pELS.Exclusive', 'others.Exclusive', 'non-cCRE.Exclusive')))

ggplot(df_to_plot, aes(x = ring, y = Frac, fill = label_2)) +
  geom_bar(stat = 'identity', width = 1, position = 'stack', size = 0.5, color = 'white') +
  coord_polar(theta = 'y') +
  theme_void() +
  theme(legend.position = 'none', axis.line.theta = element_line())





## Older code, not used
# Using GenomicRanges
cCREs_gr$OverlapCT <- FALSE
cCREs_gr$OverlapKO <- FALSE
nonPr_dREG_peaks$CT_IP$cCRE <- FALSE
nonPr_dREG_peaks$KO_IP$cCRE <- FALSE

df_to_plot <- data.frame(Genotype = c('Ctrl', 'TKO'),
                         Count = c(0, 0),
                         Fraction = c(0, 0))

idx <- findOverlaps(nonPr_dREG_peaks$CT_IP, cCREs_gr)
df_to_plot$Fraction[1] <- length(unique(queryHits(idx))) / queryLength(idx)
df_to_plot$Count[1] <- length(unique(queryHits(idx)))
cCREs_gr$OverlapCT[unique(subjectHits(idx))] <- TRUE
nonPr_dREG_peaks$CT_IP$cCRE[unique(queryHits(idx))] <- TRUE

idx <- findOverlaps(nonPr_dREG_peaks$KO_IP, cCREs_gr)
df_to_plot$Fraction[2] <- length(unique(queryHits(idx))) / queryLength(idx)
df_to_plot$Count[2] <- length(unique(queryHits(idx)))
cCREs_gr$OverlapKO[unique(subjectHits(idx))] <- TRUE
nonPr_dREG_peaks$KO_IP$cCRE[unique(queryHits(idx))] <- TRUE

nonPr_dREG_peaks$CT_IP$exclusive <- TRUE
nonPr_dREG_peaks$KO_IP$exclusive <- TRUE
idx <- findOverlaps(nonPr_dREG_peaks$CT_IP, nonPr_dREG_peaks$KO_IP)
nonPr_dREG_peaks$CT_IP$exclusive[unique(queryHits(idx))] <- FALSE
nonPr_dREG_peaks$KO_IP$exclusive[unique(subjectHits(idx))] <- FALSE

df_to_plot$Total <- c(length(nonPr_dREG_peaks$CT_IP), length(nonPr_dREG_peaks$KO_IP))

ggplot(df_to_plot, aes(x = Genotype, y = Fraction)) +
  geom_bar(stat = 'identity', fill = 'darkgreen', color = 'black', linewidth = 1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = NULL, y = 'Fraction overlapping cCREs')

KO_dREG_summary <- as.data.frame(nonPr_dREG_peaks$KO_IP)
KO_dREG_summary <- KO_dREG_summary %>% group_by(exclusive, cCRE) %>% summarise(Count = n()) %>%
  mutate(Pct = Count / sum(Count)) %>% ungroup()


cCREs_tb <- as.data.frame(cCREs_gr) %>% mutate(label = sub(",.*", "", label))
cCREs_tb <- cCREs_tb %>% subset(OverlapCT | OverlapKO) %>%
  mutate(Group = if_else(OverlapCT & OverlapKO, "Common", if_else(OverlapCT, "CT_only", "KO_only"))) %>%
  mutate(label = if_else(label != 'pELS' & label != 'dELS', 'Others', label))

df_to_plot_2 <- cCREs_tb %>% group_by(Group, label) %>% summarise(Count = n()) %>%
  mutate(Pct = Count / sum(Count)) %>% ungroup()

ggplot(df_to_plot, aes(x = Group, y = Pct, fill = label)) +
  geom_bar(stat = 'identity')


