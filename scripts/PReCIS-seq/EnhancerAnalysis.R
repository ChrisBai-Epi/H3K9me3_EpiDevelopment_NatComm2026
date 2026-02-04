library(tidyverse)
library(openxlsx)
library(ggrastr)
library(GenomicRanges)

####Read in the bed files containing dREG results
dREG_colNames <- c("chrom", "start", "end", "score", "pvalue", "center")

dRegPeaks <- list()
bed <- read.table("../dREG/CT_IP.dREG.peak.full.bed.gz", col.names = dREG_colNames)
dRegPeaks[['CT_IP']] <- GRanges(seqnames = bed$chrom,
                           ranges = IRanges(start = bed$start + 1, end = bed$end),
                           strand = rep('*', nrow(bed)),
                           score = bed$score, pvalue = bed$pvalue, center = bed$center)
bed <- read.table("../dREG/KO_IP.dREG.peak.full.bed.gz", col.names = dREG_colNames)
dRegPeaks[['KO_IP']] <- GRanges(seqnames = bed$chrom,
                           ranges = IRanges(start = bed$start + 1, end = bed$end),
                           strand = rep('*', nrow(bed)),
                           score = bed$score, pvalue = bed$pvalue, center = bed$center)
rm(bed)

dRegPeaks <- lapply(dRegPeaks, function(gr){
  seqlevels(gr) <- c(paste0("chr", 1:19), "chrX", "chrY")
  gr$name <- paste0("Peak", seq_along(gr))
  return(gr)
})


####Read in annotation and generate promoter regions, which typicall contains dREG identified peaks
Annotations <- read.table("RefSeq_MostActiveTranscriptPerGene_110123.bed", 
                          header=FALSE, stringsAsFactors = FALSE)
colnames(Annotations) <- c("chrom", "start", "end", "width", "strand", "transcript_id", "gene_name", "gene_type", "gene_id")
AllGenes <- sort(GRanges(Annotations))
CodingGenes <- subset(AllGenes, gene_type == 'protein_coding')

##Make annotation from ensembl
#AllGenes <- genes(txdb_ensembl)
#names(AllGenes) <- NULL
#seqlevelsStyle(AllGenes) <- 'UCSC'

###Remove all peaks that overlap promoter regions and also scores higher than 0.4
nonPr_dREG_peaks <- lapply(dRegPeaks, function(gr){
  pr <- promoters(AllGenes, upstream = 500, downstream = 500)
  overlapped <- findOverlaps(gr, pr)
  np_gr <- gr[setdiff(seq_along(gr), queryHits(overlapped))]
  np_gr_high <- subset(np_gr, score > 0.4)
})

local({
  idx_overlapped <- findOverlaps(nonPr_dREG_peaks$CT_IP, nonPr_dREG_peaks$KO_IP, ignore.strand = TRUE)
  overlapped_CT <- queryHits(idx_overlapped)
  overlapped_KO <- subjectHits(idx_overlapped)
  Up_CT_idx <- overlapped_CT[nonPr_dREG_peaks$CT_IP$score[overlapped_CT] > nonPr_dREG_peaks$KO_IP$score[overlapped_KO] + 0.3]
  Up_KO_idx <- overlapped_KO[nonPr_dREG_peaks$KO_IP$score[overlapped_KO] > nonPr_dREG_peaks$CT_IP$score[overlapped_CT] + 0.3]
  
  CT_exclusive_dRegPeaks <<- nonPr_dREG_peaks$CT_IP[setdiff(seq_along(nonPr_dREG_peaks$CT_IP), unique(overlapped_CT))]
  KO_exclusive_dRegPeaks <<- nonPr_dREG_peaks$KO_IP[setdiff(seq_along(nonPr_dREG_peaks$KO_IP), unique(overlapped_KO))]
  CT_Up_dRegPeaks <<- nonPr_dREG_peaks$CT_IP[Up_CT_idx]
  KO_Up_dRegPeaks <<- nonPr_dREG_peaks$KO_IP[Up_KO_idx]
})


CB_export_dREG <- function(gr, fileName, width = NULL) {
  op = options("scipen")
  options(scipen = 999)
  if (!is.null(width)){
    gr = resize(gr, width = width,fix = "center")
  }
  gr <- gr[order(as.character(seqnames(gr)), start(gr))]
  df <- as.data.frame(gr)
  df$start <- df$start - 1
  bed_df <- df[, c('seqnames', 'start', 'end', 'name', 'score', 'strand', 'pvalue', 'center')]
  write.table(bed_df, file = fileName, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  options(op)
}

CB_export_dREG(CT_exclusive_dRegPeaks, "../dREG/CT_exclusive_dREGs.bed")
CB_export_dREG(KO_exclusive_dRegPeaks, "../dREG/KO_exclusive_dREGs.bed")
CB_export_dREG(KO_Up_dRegPeaks, "../dREG/KO_Up_dREGs.bed")
CB_export_dREG(CT_Up_dRegPeaks, "../dREG/CT_Up_dREGs.bed")
CB_export_dREG(nonPr_dREG_peaks$KO_IP, "../dREG/KO_IP_all_dREGs.bed")
CB_export_dREG(nonPr_dREG_peaks$CT_IP, "../dREG/CT_IP_all_dREGs.bed")


##### Plotting methylation around dREG peaks
set.seed(428)
random_CT_dREGs <- nonPr_dREG_peaks$CT_IP[sample(length(nonPr_dREG_peaks$CT_IP), 2500)]
CB_export_dREG(random_CT_dREGs, "../dREG/Random_CT_dREGs.bed")

system2("deeptools", args = "--version")
command <- paste0("computeMatrix reference-point -S ", "../dREG/Heatmap/H3K9me3_IgG_subtracted.bw",
                  " -R ", "../dREG/Random_CT_dREGs.bed", " ", "../dREG/KO_exclusive_dREGs.bed",
                  " -o ", "../dREG/Heatmap/matrix_dREGs.gz",
                  " --referencePoint center --binSize 20 --sortRegions descend -p 15 -a 5000 -b 5000",
                  " --maxThreshold 10")
system(command)

####Need to decompress the matrix file, change the maxthreshold slots to null, and then recompress
system("gunzip ../dREG/Heatmap/matrix_dREGs.gz")
## read the data matrix and filter for the top 500 lines from both Ctrl and TKO
mat <- read.table('../dREG/Heatmap/matrix_dREGs', skip = 1)
mat_new <- mat[c(1:500, 2500:3000),]
## read the metainfo (firstline) and modify
mat_info <- read_lines('../dREG/Heatmap/matrix_dREGs', n_max = 1)
mat_info
mat_info <- paste0("@{\"upstream\":[5000],\"downstream\":[5000],\"body\":[0],\"bin size\":[20],\"ref point\":[\"center\"],\"verbose\":false,\"bin avg type\":\"mean\",\"missing data as zero\":false,",
                   "\"min threshold\":null,\"max threshold\":null,\"scale\":1,\"skip zeros\":false,\"nan after end\":false,\"proc number\":15,\"sort regions\":\"descend\",\"sort using\":\"mean\",\"unscaled 5 prime\":[0],\"unscaled 3 prime\":[0],\"sample_labels\":[\"H3K9me3_IgG_subtracted\"],\"group_labels\":[\"Random_CT_dREGs.bed\",\"KO_exclusive_dREGs.bed\"],",
                   "\"sample_boundaries\":[0,500],\"group_boundaries\":[0,500,1001]}")

## Write information into a new matrix file for plotting
write_lines(mat_info, file = '../dREG/Heatmap/matrix_dREGs_filtered')
write.table(mat_new, file = '../dREG/Heatmap/matrix_dREGs_filtered', row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = '\t', append = TRUE)
rm(mat, mat_info, mat_new)

system("gzip ../dREG/Heatmap/matrix_dREGs_filtered")

command <- paste0("plotHeatmap -m ", "../dREG/Heatmap/matrix_dREGs_filtered.gz",
                 " -o ", "../dREG/Heatmap/Methylation_dREGs.pdf",
                 " --colorMap RdYlBu_r --samplesLabel H3K9me3 --regionsLabel Ctrl TKO --refPointLabel PeakCenter",
                 " --heatmapWidth 8 --heatmapHeight 16 --averageTypeSummaryPlot mean --zMin -0.2 --zMax 0.6")
system(command)

system("rm ../dREG/Heatmap/matrix_dREGs*")



######## Prepare files for homer motif analysis ########
CB_PrepareBed_forHomer <- function(gr, fileName) {
  gr <- gr[order(seqnames(gr), start(gr))]
  df <- as.data.frame(gr)
  df$start <- df$center - 150
  df$end <- df$center + 150
  df$name <- paste0("Peak", row.names(df))
  bed_df <- df[, c('seqnames', 'start', 'end', 'name', 'score', 'strand')]
  write.table(bed_df, file = fileName, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

CB_PrepareBed_forHomer(nonPr_dREG_peaks$CT_IP, "../dREG/ForHomer_dREG_CT_all.bed")
CB_PrepareBed_forHomer(nonPr_dREG_peaks$KO_IP, "../dREG/ForHomer_dREG_KO_all.bed")
CB_PrepareBed_forHomer(CT_exclusive_dRegPeaks, "../dREG/ForHomer_dREG_CT_exclusive.bed")
CB_PrepareBed_forHomer(KO_exclusive_dRegPeaks, "../dREG/ForHomer_dREG_KO_exclusive.bed")

set.seed(675)
random_CT_dREGs <- nonPr_dREG_peaks$CT_IP[sample(length(nonPr_dREG_peaks$CT_IP), 10000)]
CB_PrepareBed_forHomer(random_CT_dREGs, "../dREG/ForHomer_dREG_Random_CT_15k.bed")

rm(CB_PrepareBed_forHomer)


####### Associate enhancer to genes #########################
AllGenes$Width <- width(AllGenes)
TSS_regions <- promoters(AllGenes, upstream = 5, downstream = 5)
TSS_regions <- TSS_regions[order(seqnames(TSS_regions), start(TSS_regions))]

local({
  df <- as.data.frame(TSS_regions)
  df$start <- df$start - 1
  df$score <- 0
  bed_df <- df[, c('seqnames', 'start', 'end', 'gene_id', 'Width', 'strand', 'gene_name', 'gene_type')]
  write.table(bed_df, file = '../dREG/TSS_codingGenes.bed', sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
})

system("bedtools --version")
command <- paste0("bedtools closest -a ../dREG/TSS_codingGenes.bed ",
                  "-b ../dREG/KO_exclusive_dREGs.bed ",
                  "-D a -k 5 ",
                  "1> ../dREG/dREG_KO_exclusive_distribution.bed 2> /dev/null")
system(command)

command <- paste0("bedtools closest -a ../dREG/TSS_codingGenes.bed ",
                  "-b ../dREG/KO_IP_all_dREGs.bed ",
                  "-D a -k 5 ",
                  "1> ../dREG/dREG_KO_all_distribution.bed 2> /dev/null")
system(command)


###### Analyze obtained enhancer-promoter distance table #######

PreciseSeqFC_table <- read.xlsx("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/Precise2024_FoldChangeTable_annotated.xlsx")
PreciseSeqFC_table$Category <- factor(PreciseSeqFC_table$Category, levels = c('upPI_nsGB', 'nsPI_upGB', 'downPI_upGB', 'downPI_nsGB', 'others'))
Nascent_GOI <- list()
for (name in unique(PreciseSeqFC_table$Category)){
  Nascent_GOI[[name]] <- subset(PreciseSeqFC_table, subset = Category == name, select = GeneID, drop = TRUE)
  write.table(Nascent_GOI[[name]], file = paste0('MotifAnalysisFiles/PreciseSeq_', name, '.txt'), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

Unregulated_genes <- subset(PreciseSeqFC_table, subset = Reg_PI == 'ns' & Reg_GB == 'ns', select = GeneID, drop = TRUE)
set.seed(428)
Nascent_GOI[['Others']] <- sample(Unregulated_genes, size = 200)

CB_Prepare_dREG_Distribution_df <- function(fileName){
  
  dREG_Dist <- read.table(fileName) %>% select(!c("V14", "V16"))
  colnames(dREG_Dist) <- c("chrom", "start", "end", "ensembl", "Width", "strand", "name", "biotype", 
                           "peak_chrom", "peak_start", "peak_end", "peak_name", "score", "pvalue", "distance")
  
  dREG_Dist <- dREG_Dist %>%
    mutate(category = case_when(ensembl %in% Nascent_GOI$Up_PI_ns_GB_genes ~ "upPI_nsGB",
                                ensembl %in% Nascent_GOI$Down_PI_ns_GB_genes ~ "downPI_nsGB",
                                ensembl %in% Nascent_GOI$Ns_PI_Up_GB_genes ~ "nsPI_upGB",
                                ensembl %in% Nascent_GOI$Down_PI_Up_GB_genes ~ "downPI_upGB",
                                .default = "Others")) %>%
    mutate(position = if_else(abs(distance) < Width & distance < 0, "Intronic", "Intergenic")) %>%
    mutate(distance = if_else(abs(distance) > Width, abs(distance), distance))
  dREG_df <- dREG_Dist %>% select(!c(peak_chrom:pvalue)) %>%
    mutate(logDist = case_when(distance < 0 ~ -log10(-distance),
                               distance == 0 ~ 0,
                               distance > 0 ~ log10(distance))) %>%
    group_by(ensembl) %>% arrange(distance, .by_group = TRUE) %>% 
    mutate(Dist_rank = row_number()) %>%
    mutate(closest = min(logDist)) %>% ungroup()
  dREG_df$Dist_rank <- factor(dREG_df$Dist_rank)
  dREG_df$category <- factor(dREG_df$category, levels = c("upPI_nsGB", "nsPI_upGB", "downPI_upGB", "downPI_nsGB", "Others"))
  return(dREG_df)
}

CB_plot_dREG_Number <- function(df, range){
  dREG_numbers_forGene <- df %>% group_by(ensembl, name, category) %>%
    summarise(Num = sum(abs(distance) < range))
  dREG_numbers_forCategory <- dREG_numbers_forGene %>% group_by(category) %>%
    summarise("1" = mean(Num == 1), "2" = mean(Num == 2), "3" = mean(Num == 3),
              "4" = mean(Num == 4), ">=5" = mean(Num == 5)) %>%
    pivot_longer(cols = 2:6, names_to = "Num_dREGs", values_to = "Fraction")
  dREG_numbers_forCategory$Num_dREGs <- factor(dREG_numbers_forCategory$Num_dREGs, levels = c('>=5', '4', '3', '2', '1'))
  
  ggplot(dREG_numbers_forCategory, aes(x = category, y = Fraction, fill = Num_dREGs)) +
    geom_bar(stat = 'identity')
}

df_to_plot <- CB_Prepare_dREG_Distribution_df("../dREG/dREG_CT_all_distribution.bed")
df_to_plot <- CB_Prepare_dREG_Distribution_df("../dREG/dREG_KO_all_distribution.bed")
df_to_plot <- CB_Prepare_dREG_Distribution_df("../dREG/dREG_KO_exclusive_distribution.bed")

ggplot(subset(df_to_plot, Dist_rank == 1), aes(x = closest, colour = category)) +
  stat_ecdf() +
  labs(x = "Log10(Distance to TSS)",
       y = "Cumulative Percentage")

ggplot(df_to_plot, aes(x = category, y = logDist, colour = Dist_rank)) +
  geom_boxplot() +
  facet_wrap(~position) +
  labs(x = element_blank(),
       y = "Log10(Distance to TSS)")

CB_plot_dREG_Number(df_to_plot, 1e4) + 
  labs(x = element_blank(), fill = "Num of dREGs within 10k") +
  theme(legend.title = element_text(angle = 90), legend.title.position = 'left')

##### Analyze genes associated with enhancers
CB_plot_dREG_associated_Genes <- function(df, range){
  dREG_numbers_forGene <- df %>% group_by(ensembl, name, category) %>%
    summarise(Num = sum(abs(distance) < range))
  Gene_Category_for_dREGs <- dREG_numbers_forGene %>% group_by(Num, category) %>%
    summarise(Count = n()) %>% mutate(Percent = Count / sum(Count)) %>% ungroup() %>%
    subset(category != 'Others')
  
  Gene_Category_for_dREGs$Num <- factor(Gene_Category_for_dREGs$Num, levels = c('0', '1', '2', '3', '4', '5'))
  
  ggplot(Gene_Category_for_dREGs, aes(x = Num, y = Percent, fill = category)) +
    geom_bar(stat = 'identity') +
    scale_y_continuous(labels = scales::percent)
}

df_to_plot <- CB_Prepare_dREG_Distribution_df("../dREG/dREG_KO_exclusive_distribution.bed")
CB_plot_dREG_associated_Genes(df_to_plot, 1e4) +
  labs(x = 'Num of dREG peaks within 10kb')



################
## Only analyzing the closest dReg peak, similar to Epic2 H3K9me3 peak analysis
## Prepare for Gene-closest-dReg file
TSS_regions <- promoters(AllGenes, upstream = 500, downstream = 500, use.names = TRUE)
TSS_regions <- subset(TSS_regions, start(TSS_regions) >= 0)
TSS_regions <- TSS_regions[order(seqnames(TSS_regions), start(TSS_regions))]

##Check KO_exclusive_dREGs and export as bed files; obtaining process see Enhancer_plotting_dREGs.R
KO_exclusive_dREGs

local({
  gr <- KO_exclusive_dREGs
  gr <- gr[order(seqnames(gr), start(gr))]
  df <- as.data.frame(gr)
  df$start <- df$start - 1
  df$score <- df$methyl_signal
  bed_df <- df[, c('seqnames', 'start', 'end', 'name', 'score', 'strand', 'MethylTag')]
  write.table(bed_df, file = '../dREG/dREG_KO_exclusive_withMethylInfo', sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
})

local({
  df <- as.data.frame(TSS_regions)
  df$start <- df$start - 1
  df$score <- 0
  bed_df <- df[, c('seqnames', 'start', 'end', 'gene_id', 'score', 'strand', 'gene_name', 'gene_type')]
  write.table(bed_df, file = '../dREG/TSS_500bp_allGenes.bed', sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
})

system("bedtools --version")
command <- paste0("bedtools closest -a ", "../dREG/TSS_500bp_allGenes.bed",
                  " -b ", "../dREG/ForHomer_dREG_CT_all.bed",
                  " -D a -t first > ",
                  "../dREG/Closest_CTall_dRegPeaktoPromoters.bed")
system(command)
command <- paste0("bedtools closest -a ", "../dREG/TSS_500bp_allGenes.bed",
                  " -b ", "../dREG/ForHomer_dREG_KO_all.bed",
                  " -D a -t first > ",
                  "../dREG/Closest_KOall_dRegPeaktoPromoters.bed")
system(command)
command <- paste0("bedtools closest -a ", "../dREG/TSS_500bp_allGenes.bed",
                  " -b ", "../dREG/ForHomer_dREG_KO_exclusive.bed",
                  " -D a -t first > ",
                  "../dREG/Closest_KOexclu_dRegPeaktoPromoters.bed")
system(command)
command <- paste0("bedtools closest -a ", "../dREG/TSS_500bp_allGenes.bed",
                  " -b ", "../dREG/dREG_KO_exclusive_withMethylInfo",
                  " -D a -t first > ",
                  "../dREG/Closest_KOexclu_dRegToTSS_withMethyl.bed")
system(command)

## For dReg_distance, pick either file for KO_all or KO_exclusive
dReg_distance <- read.table("../dREG/Closest_KOexclu_dRegPeaktoPromoters.bed", sep = "\t") %>%
  select(c(V4, V7, V8, V13, V15)) %>%
  rename(GeneID = V4, Symbol = V7, Biotype = V8, PeakScore = V13, Distance = V15) %>%
  filter(PeakScore > 0)

## Tag upregulated genes based on either bulk RNA-seq data or nascent DE data
df_to_plot <- dReg_distance %>% mutate(Group = case_when(GeneID %in% RNAseq_UpGenes ~ "Upregulated",
                                                            .default =  "Others")) %>%
  mutate(Group = factor(Group, levels = c('Upregulated', 'Others'))) %>% arrange(desc(Group))

df_to_plot <- dReg_distance %>% mutate(Group = case_when(GeneID %in% Nascent_UpGenes ~ "Upregulated",
                                                         .default =  "Others")) %>%
  mutate(Group = factor(Group, levels = c('Upregulated', 'Others'))) %>% arrange(desc(Group))



ggplot(df_to_plot, aes(x = Group, y = abs(Distance))) +
  geom_bar(stat = 'summary', fun = 'mean', fill = 'darkgreen', width = 0.8, color = 'black', size = 1) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, size =1) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = element_blank(), y = "Mean distance to closest dREG peaks") +
  theme(axis.text.x = element_text(size = 12))

ggplot(df_to_plot, aes(x = Group, y = abs(Distance))) +
  geom_boxplot(outliers = FALSE, fill = 'darkgreen', color = 'black', linewidth = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = element_blank(), y = "Mean distance to closest dREG peaks") +
  theme(axis.text.x = element_text(size = 12))

#Perform t test for above plot
a <- subset(df_to_plot, Group == 'Upregulated')
b <- subset(df_to_plot, Group == 'Others')
t.test(abs(a$Distance), abs(b$Distance))


###### Plot log2FC vs distance ######
## Use either bulk RNA-seq data or Nascent DE data
RNAseq_E16_DEs <- read.xlsx('/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/RNAseq_E16_genotypeOnly_noH614_DEgenes_annotated.xlsx', 
                            sheet = 1) %>% distinct(ensembl, .keep_all = TRUE) %>%
  select(c(ensembl, log2FoldChange, padj)) %>% rename(GeneID = ensembl)

Nascent_fullGene_DE <- read.xlsx("DiffExpressedGenes_fullGene.xlsx", rowNames = TRUE) %>%
  distinct(Ensembl, .keep_all = TRUE) %>%
  select(c(Ensembl, log2FoldChange, padj)) %>% rename(GeneID = Ensembl)

dReg_distance <- read.table("../dREG/Closest_KOexclu_dRegToTSS_withMethyl.bed", sep = "\t") %>%
  select(c(V4, V7, V8, V13, V15, V16)) %>%
  rename(GeneID = V4, Symbol = V7, Biotype = V8, MethylSignal = V13, Tag = V15, Distance = V16) %>%
  filter(Tag != '.') %>%
  inner_join(Nascent_fullGene_DE, by = 'GeneID')

df_to_plot <- dReg_distance %>%
  mutate(Group = if_else(padj < 0.05 & log2FoldChange > 1, 'Upregulated', 'Others')) %>%
  mutate(Group = factor(Group, levels = c('Upregulated', 'Others'))) %>% arrange(desc(Group)) %>%
  mutate(LogDistance = case_when(Distance > 0 ~ log10(Distance),
                                 Distance == 0 ~ 0,
                                 Distance < 0 ~ -log10(-Distance)))

ggplot(df_to_plot, aes(x = Distance, y = log2FoldChange, color = Group)) +
  rasterize(geom_point(size = 1), dpi = 300) +
  scale_color_manual(values = c('red', 'gray')) +
  scale_x_continuous(limits = c(-1e6, 1e6)) +
  labs(x = 'Distance from TSS to closest dREG', y = 'Log2FC(TKO/Ctrl)', color = NULL)

ggplot(df_to_plot, aes(x = abs(LogDistance), color = Group)) +
  stat_ecdf() +
  labs(x = 'Log10(Distance)', y = 'Cumulative Percentage')

# Find distribution of TSS-to-dREG distances
distance_distribution <- data.frame(Bins = c('0-10kb', '10-20kb', '20-50kb', '50-200kb', '200-500kb', '>500kb'))
genes <- subset(df_to_plot, Group == 'Upregulated')
genes <- df_to_plot
Bins <- cut(abs(genes$Distance), breaks = c(0, 10000, 20000, 50000, 200000, 500000, Inf),
            labels = c('0-10kb', '10-20kb', '20-50kb', '50-200kb', '200-500kb', '>500kb'))
Bins <- table(Bins)
distance_distribution$Upregulated <- Bins / sum(Bins)
distance_distribution$allGenes <- Bins / sum(Bins)
distance_distribution <- pivot_longer(distance_distribution, cols = c(Upregulated, allGenes), names_to = 'Group', values_to = 'Frac')
ggplot(distance_distribution, aes(x = Bins, y = Frac, fill = Group)) +
  geom_bar(stat = 'identity', position = position_dodge())

# Plot with PreciseSeq gene lists
PreciseSeqFC_table <- read.xlsx("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/Precise2024_FoldChangeTable_annotated.xlsx")
PreciseSeqFC_table$Category <- factor(PreciseSeqFC_table$Category,
                                      levels = c('upPI_nsGB', 'nsPI_upGB', 'downPI_upGB', 'downPI_nsGB', 'others'))

Nascent_GOI <- list()
for (name in unique(PreciseSeqFC_table$Category)){
  Nascent_GOI[[name]] <- subset(PreciseSeqFC_table, subset = Category == name, select = GeneID, drop = TRUE)
  write.table(Nascent_GOI[[name]], file = paste0('MotifAnalysisFiles/PreciseSeq_', name, '.txt'), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

dReg_distance <- read.table("../dREG/Closest_KOexclu_dRegPeaktoPromoters.bed", sep = "\t") %>%
  select(c(V4, V7, V8, V13, V15)) %>%
  rename(GeneID = V4, Symbol = V7, Biotype = V8, PeakScore = V13, Distance = V15) %>%
  filter(PeakScore > 0) %>%
  inner_join(PreciseSeqFC_table[,c('GeneID', 'DEseq_log2FC_PI', 'DEseq_log2FC_GB', 'Category')], by = 'GeneID')

df_to_plot <- dReg_distance %>% arrange(desc(Category))

ggplot(df_to_plot, aes(x = Distance, y = DEseq_log2FC_GB, color = Category)) +
  rasterize(geom_point(size = 1), dpi = 300) +
  scale_color_manual(values = c('red', 'black', 'blue', 'purple', 'gray')) +
  scale_x_continuous(limits = c(-1e6, 1e6)) +
  labs(x = 'Distance from TSS to closest dREG', y = 'GB Log2FC(TKO/Ctrl)', color = NULL)

ggplot(df_to_plot, aes(x = Distance, y = DEseq_log2FC_PI, color = Category)) +
  rasterize(geom_point(size = 1), dpi = 300) +
  scale_color_manual(values = c('red', 'black', 'blue', 'purple', 'gray')) +
  scale_x_continuous(limits = c(-1e6, 1e6)) +
  labs(x = 'Distance from TSS to closest dREG', y = 'PI Log2FC(TKO/Ctrl)', color = NULL)

ggplot(df_to_plot, aes(x = Category, y = abs(Distance))) +
  geom_bar(stat = 'summary', fun = 'mean', fill = 'darkgreen', width = 0.8, color = 'black', size = 1) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, size =1) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = element_blank(), y = "Mean distance to closest H3K9me3 peaks") +
  theme(axis.text.x = element_text(size = 12))
