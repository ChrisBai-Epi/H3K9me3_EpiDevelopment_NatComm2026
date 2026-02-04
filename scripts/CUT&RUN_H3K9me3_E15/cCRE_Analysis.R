library(tidyverse)
library(GenomicRanges)

#Load the peak file
Epic2_peaks <- read.table("DataFiles/epic2_peaks_H3K9me3_E15.bed", header = TRUE)
Epic2_peaks$PeakID <- paste0("Epic", row.names(Epic2_peaks))
Epic2_peaks <- Epic2_peaks[, c(1, 2, 3, 11, 5, 6, 7, 8, 9, 10, 4)]

write.table(Epic2_peaks[, 1:6], "DataFiles/Epic2_peaks.bed", col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
Epic2_peaks <- GRanges(Epic2_peaks)

cCREs <- read.table("DataFiles/mm10_cCREs_all_ENCFF904ZZH.bed")
colnames(cCREs) <- c("seqnames", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRGB", "Source", "Type")
cCREs <- GRanges(cCREs)


#Find out the number of methylated cCRE region around every gene
local({
  idx <- findOverlaps(query = Epic2_peaks, subject = cCREs)
  methylated_cCREs <<- cCREs[unique(subjectHits(idx))]
})

local({
  pr <- promoters(all_genes, upstream = 10000, downstream = 10000)
  num_overlaps <- countOverlaps(pr, methylated_cCREs)
  all_genes$Num_pr_methy_cCREs <<- num_overlaps
})

local({
  expanded_gr <- resize(all_genes, width = width(all_genes) + 20000, fix = 'center')
  num_overlaps <- countOverlaps(expanded_gr, methylated_cCREs)
  all_genes$Num_flk_methy_cCREs <<- num_overlaps
})

Num_cCRE_allgenes <- as.data.frame(all_genes)
Num_cCRE_GOI <- lapply(list_GOI, function(lt){
  df <- Num_cCRE_allgenes %>% subset(gene_id %in% lt | Symbol %in% lt)
})

df_to_plot <- do.call(rbind, Num_cCRE_GOI) %>% rownames_to_column(var = "GeneSet")
df_to_plot$GeneSet <- factor(sub("\\..*", "", df_to_plot$GeneSet), levels = c("RNAseq_up", "RNAseq_down", "SC_upGenes", "SC_downGenes", "Precise_up", "Random"))

df_to_plot <- df_to_plot %>% group_by(GeneSet) %>%
  summarise(mean_num_pr = mean(Num_pr_methy_cCREs), frac_pr_w_cCREs = mean(Num_pr_methy_cCREs > 0),
            mean_num_flk = mean(Num_flk_methy_cCREs), frac_flk_w_cCREs = mean(Num_flk_methy_cCREs > 0)) %>%
  ungroup() %>% subset(GeneSet != 'Precise_up')

ggplot(df_to_plot, aes(x = GeneSet, y = mean_num_pr)) +
  geom_bar(stat = 'identity', fill = 'forestgreen') +
  labs(x = element_blank(), y = "Mean number of methylated cCREs")

ggplot(df_to_plot, aes(x = GeneSet, y = mean_num_flk)) +
  geom_bar(stat = 'identity', fill = 'forestgreen') +
  labs(x = element_blank(), y = "Mean number of methylated cCREs")







