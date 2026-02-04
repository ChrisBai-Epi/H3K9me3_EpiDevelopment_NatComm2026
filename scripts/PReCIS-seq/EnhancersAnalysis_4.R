library(tidyverse)
library(openxlsx)
library(GenomicRanges)

## This script focus on finding out whether methylation around genes contains active dREGs
## Also checks the analyze distance between dREGs and activated genes

## Analyze nearest Epic2 peak (from TSS) for each gene
Epic2_peaks <- read.table("../BrowserTracks/Epic2_peaks.bed")
colnames(Epic2_peaks) <- c('seqnames', 'start', 'end', 'name', 'score', 'strand')
Epic2_peaks <- GRanges(Epic2_peaks)

Peak_distance <- read.table('/Users/kb744/Tumbar Lab Local/CUTandRUN/Functional_Analysis/DataFiles/Methyl_Annotation/ClosestEpic2Peak_toPromoters.bed', sep = "\t") %>%
  select(c(V4, V10, V11, V13)) %>% rename(ensembl = V4, PeakScore = V11, Distance = V13, PeakID = V10) %>% filter(PeakScore > 0)


Epic2_peaks$dREG_CT <- FALSE
idx <- findOverlaps(Epic2_peaks, nonPr_dREG_peaks$CT_IP)
Epic2_peaks$dREG_CT[unique(queryHits(idx))] <- TRUE

Epic2_peaks$dREG_KO <- FALSE
idx <- findOverlaps(Epic2_peaks, nonPr_dREG_peaks$KO_IP)
Epic2_peaks$dREG_KO[unique(queryHits(idx))] <- TRUE

Epic2_peaks$cCRE <- FALSE
idx <- findOverlaps(Epic2_peaks, cCREs_gr, select = 'first')
Epic2_peaks$cCRE <- cCREs_gr$label[idx]

cCREs_gr$Epic2 <- FALSE
idx <- findOverlaps(Epic2_peaks, cCREs_gr)
cCREs_gr$Epic2[unique(subjectHits(idx))] <- TRUE


## Define gene lists of interest
Nascent_up_genes <- subset(Nascent_fullGene_DE, log2FoldChange > 1 & padj < 0.05, select = 'Ensembl', drop = TRUE)
set.seed(428)
Nascent_random_genes <- sample(Nascent_fullGene_DE$Ensembl, 500)
Nascent_regulated_genes <- subset(PreciseSeqFC_table, Category != 'others', select = 'GeneID', drop = TRUE)

## Check the relationship between different gene lists
library(ggvenn)
ggvenn(list("Nascent_up" = Nascent_up_genes, "Nascent_regulated" = Nascent_regulated_genes))


## Define a function to subset for epic2 peaks within a range around genes of interest
all_promoters <- promoters(AllGenes, upstream = 50000, downstream = 50000)

chrs <- c(paste0('chr', 1:19), 'chrX', 'chrY', 'chrM')
all_promoters <- keepSeqlevels(all_promoters, chrs, pruning.mode = 'tidy')
Epic2_peaks <- keepSeqlevels(Epic2_peaks, chrs[-22], pruning.mode = 'tidy')

CB_get_peaks_around_genes <- function(peaks, gList){
  require(GenomicRanges)
  
  gr <- subset(all_promoters, gene_id %in% gList)
  idx <- findOverlaps(gr, peaks)
  peak_of_interest <- peaks[unique(subjectHits(idx))]
  
}

epic_nascent_up_genes <- CB_get_peaks_around_genes(Epic2_peaks, Nascent_up_genes)
epic_nascent_random_genes <- CB_get_peaks_around_genes(Epic2_peaks, Nascent_random_genes)

sum(epic_nascent_up_genes$dREG_KO)
sum(epic_nascent_random_genes$dREG_KO)

rm(CB_get_peaks_around_genes, all_promoters)

## Plot types of cCREs overlapping Epic2 peaks (methylated)
Methylated_cCREs <- subset(cCREs_gr, Epic2)
Methy_cCRE_near_genes <- CB_get_peaks_around_genes(Methylated_cCREs, all_promoters$gene_id)
Methylated_cCREs$nearGene <- FALSE
Methylated_cCREs$nearGene[match(Methy_cCRE_near_genes$name, Methylated_cCREs$name)] <- TRUE

df_to_plot <- as.data.frame(Methylated_cCREs) %>% group_by(label) %>% summarise(Count = n()) %>%
  mutate(Pct = Count / sum(Count))

ggplot(df_to_plot, aes(x = '', y = Pct, fill = label)) +
  geom_bar(stat = 'identity') +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_fill_brewer(palette = 'Set3') +
  labs(x = NULL, y = "Percentage")

### Check genes that have epic2-covered dREG peaks around their promoter regions
epic_regulatory <- subset(Epic2_peaks, dREG_KO)

idx <- findOverlaps(all_promoters, epic_regulatory)
H3K9_regulated_genes <- all_promoters[unique(queryHits(idx))]
set.seed(428)
random_genes <- all_promoters[sample(seq_along(all_promoters), 500)]

length(intersect(Nascent_up_genes, H3K9_regulated_genes$gene_id))
length(intersect(Nascent_up_genes, random_genes$gene_id))

rm(list = ls()[grep("^epic_", ls())])
rm(list = ls()[grep("^H3K9_", ls())])
rm(all_promoters)



