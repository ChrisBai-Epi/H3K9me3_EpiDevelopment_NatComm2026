library(tidyverse)
library(openxlsx)

#### Get annatation for all genes ####
library(GenomicFeatures)
library(GenomicRanges)
library(GenomeInfoDb)
library(txdbmaker)

txdb_ensembl <- makeTxDbFromEnsembl(organism = 'Mus musculus', release = 110)
all_genes <- genes(txdb_ensembl)
info <- select(Ens_mm39_110, keys = unlist(all_genes$gene_id), keytype = 'GENEID', 
               columns = c('GENEBIOTYPE', 'GENEID', 'GENENAME', 'SYMBOL', 'DESCRIPTION'))
all_genes$Biotype <- info$GENEBIOTYPE
all_genes$Symbol <- info$SYMBOL
rm(info, txdb_ensembl)

chrs <- c(paste0('chr', 1:19), 'chrX', 'chrY', 'chrM')
seqlevelsStyle(all_genes) <- 'UCSC'
all_genes <- keepSeqlevels(all_genes, chrs, pruning.mode = 'tidy')
all_genes <- sort(all_genes)
names(all_genes) <- NULL

#### Identify genes near Epic2 peaks ####
all_promoters <- promoters(all_genes, upstream = 500, downstream = 500, use.names = TRUE)
all_promoters <- subset(all_promoters, start(all_promoters) >= 0)

all_genes_extended_10kb <- resize(all_genes, width = width(all_genes) + 20000, fix = 'center')

FoldChangeTable_Epic2_peaks <- read.xlsx('DataFiles/FoldChangeTable_Epic2_peaks_v1.xlsx') %>%
  mutate(Dynamics = factor(Dynamics, levels = c('Constant', 'Dynamic_1', 'Dynamic_2', 'Dynamic_3', 'Dynamic_4')))

identical(Epic2_union_peak$PeakID, FoldChangeTable_Epic2_peaks$PeakID)
Epic2_union_peak$Dynamics <- FoldChangeTable_Epic2_peaks$Dynamics

GR_epic2_union_peaks <- keepSeqlevels(GRanges(Epic2_union_peak), chrs[1:21], pruning.mode = 'tidy')
GR_epic2_union_peaks$Dynamics <- factor(GR_epic2_union_peaks$Dynamics, levels = c('Constant', 'Dynamic_1', 'Dynamic_2', 'Dynamic_3', 'Dynamic_4'))

Methyl_associated_genes <- list()
for (name in levels(GR_epic2_union_peaks$Dynamics)){
  cnt <- countOverlaps(all_genes_extended_10kb, subset(GR_epic2_union_peaks, Dynamics == name))
  peaks <- all_genes_extended_10kb[cnt > 0, ]
  Methyl_associated_genes[[name]] <- peaks$gene_id
}
rm(name, cnt, peaks)
rm(all_genes_extended_10kb)

#### Perform GO analysis ####
library(clusterProfiler)
library(org.Mm.eg.db)

GO_results <- lapply(Methyl_associated_genes, function(glist){
  go <- enrichGO(glist, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "ENSEMBL", pvalueCutoff = 0.05, readable = TRUE)
})

a <- lapply(GO_results, function(go){
  df <- as.data.frame(go)
})

barplot(GO_results[['Constant']], showCategory = 15)

wb <- createWorkbook()
for (name in names(GO_results)){
  cat(name, '\n')
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = as.data.frame(GO_results[[name]]))
}
saveWorkbook(wb, file = 'DataFiles/GO_DynamicPeaks.xlsx')
rm(wb, name)


#### Plot GO for Figures ####
library(patchwork)

#### Dynamic 2: increasing peaks
go_table <- as.data.frame(GO_results[['Dynamic_2']])
GO_terms_to_plot <- c("GO:0021700", "GO:0030098", "GO:0045165", "GO:0016358", "GO:0048663",
                      "GO:0033002", "GO:0048562", "GO:0008544", "GO:0007498", "GO:0030900")

df_to_plot <- go_table %>% arrange(Count) %>% mutate(Description = factor(Description, levels = Description)) %>%
  mutate(LogPadjust = -log10(p.adjust)) %>%
  subset(ID %in% GO_terms_to_plot)

p1 <- 
ggplot(df_to_plot, aes(x = Description, y = Count, fill = LogPadjust)) +
  geom_bar(stat = 'identity', color = 'black') +
  coord_flip() +
  scale_fill_gradient(low = 'white', high = 'darkgreen', limits = c(0, 8)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(x = NULL, y = "Gene count", fill = "-log10(p.adjust)", title = 'Dynamic 2') +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12))


#### Dynamic 3: decreasing peaks
go_table <- as.data.frame(GO_results[['Dynamic_3']])
GO_terms_to_plot <- c("GO:0050807", "GO:0051963", "GO:0008038", "GO:0098982",
                      "GO:0098742", "GO:0005912")

df_to_plot <- go_table %>% arrange(Count) %>% mutate(Description = factor(Description, levels = Description)) %>%
  mutate(LogPadjust = -log10(p.adjust)) %>%
  subset(ID %in% GO_terms_to_plot)

p2 <-
ggplot(df_to_plot, aes(x = Description, y = Count, fill = LogPadjust)) +
  geom_bar(stat = 'identity', color = 'black') +
  coord_flip() +
  scale_fill_gradient(low = 'white', high = 'darkgreen', limit = c(0, 8), na.value = 'darkgreen') +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(x = NULL, y = "Gene count", fill = "-log10(p.adjust)", title = 'Dynamic 3') +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12))

p1/p2 + plot_layout(heights = c(5, 3))

rm(go_table, GO_terms_to_plot, p1, p2)

#### Plot for dynamic 1, dynamic 4 and constant peaks
#### Dynamic 1: increasing peaks
go_table <- as.data.frame(GO_results[['Dynamic_1']])
GO_terms_to_plot <- c("GO:0061515", "GO:0099173", "GO:0031589", "GO:0030111",
                      "GO:0005604", "GO:0031012")

df_to_plot <- go_table %>% arrange(Count) %>% mutate(Description = factor(Description, levels = Description)) %>%
  mutate(LogPadjust = -log10(p.adjust)) %>%
  subset(ID %in% GO_terms_to_plot)

p1 <- 
  ggplot(df_to_plot, aes(x = Description, y = Count, fill = LogPadjust)) +
  geom_bar(stat = 'identity', color = 'black') +
  coord_flip() +
  scale_fill_gradient(low = 'white', high = 'darkgreen', limits = c(0, 8)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(x = NULL, y = "Gene count", fill = "-log10(p.adjust)", title = 'Dynamic 1') +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12))

#### Dynamic 4: decreasing peaks
go_table <- as.data.frame(GO_results[['Dynamic_4']])
GO_terms_to_plot <- c("GO:0060078", "GO:0008038", "GO:0007215", "GO:0042101")

df_to_plot <- go_table %>% arrange(Count) %>% mutate(Description = factor(Description, levels = Description)) %>%
  mutate(LogPadjust = -log10(p.adjust)) %>%
  subset(ID %in% GO_terms_to_plot)

p2 <- 
  ggplot(df_to_plot, aes(x = Description, y = Count, fill = LogPadjust)) +
  geom_bar(stat = 'identity', color = 'black') +
  coord_flip() +
  scale_fill_gradient(low = 'white', high = 'darkgreen', limits = c(0, 8)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(x = NULL, y = "Gene count", fill = "-log10(p.adjust)", title = 'Dynamic 4') +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12))

#### Constant Peaks
go_table <- as.data.frame(GO_results[['Constant']])
GO_terms_to_plot <- c("GO:0019236", "GO:0098742", "GO:0035107", "GO:0021536", "GO:0055123",
                      "GO:0048665", "GO:0099173", "GO:0048706", "GO:0050839")

df_to_plot <- go_table %>% arrange(Count) %>% mutate(Description = factor(Description, levels = Description)) %>%
  mutate(LogPadjust = -log10(p.adjust)) %>%
  subset(ID %in% GO_terms_to_plot)

p3 <- 
  ggplot(df_to_plot, aes(x = Description, y = Count, fill = LogPadjust)) +
  geom_bar(stat = 'identity', color = 'black') +
  coord_flip() +
  scale_fill_gradient(low = 'white', high = 'darkgreen', limits = c(0, 8), na.value = 'darkgreen') +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(x = NULL, y = "Gene count", fill = "-log10(p.adjust)", title = 'Constant') +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12))

p4 <- p1/p2 + plot_layout(heights = c(6, 4))
p3 + p4

rm(go_table, GO_terms_to_plot, p1, p2, p3, p4)
