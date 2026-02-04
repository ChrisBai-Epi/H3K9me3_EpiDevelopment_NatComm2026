library(openxlsx)
library(org.Mm.eg.db)

#################### Compare with E16 bulk RNA-seq data ########
PreciseSeqFC_table <- read.xlsx("/Users/kb744/Tumbar Lab Local/Precise-seq_2024/ExpressionAnalysis/FoldChangeTable_annotated.xlsx")
PreciseSeqFC_table$Category <- factor(PreciseSeqFC_table$Category, levels = c('upPI_nsGB', 'nsPI_upGB', 'downPI_upGB', 'downPI_nsGB', 'others'))

Nascent_GOI <- list()
for (name in unique(PreciseSeqFC_table$Category)){
  Nascent_GOI[[name]] <- subset(PreciseSeqFC_table, subset = Category == name, select = GeneID, drop = TRUE)
}

df_to_plot <- merge(PreciseSeqFC_table, RNAseq_E16_DEs, by = 'GeneID', all.x = TRUE) %>%
  rename(RNAseq_padj = padj.y, RNAseq_log2FC = log2FoldChange)

ggplot(df_to_plot, aes(x = DEseq_log2FC_GB, y = RNAseq_log2FC)) +
  geom_point() +
  geom_hline(yintercept = 1, color = 'blue', linetype = 'dashed') +
  geom_hline(yintercept = 0, color = 'black', linetype = 'solid') +
  geom_vline(xintercept = 1, color = 'blue', linetype = 'dashed') +
  geom_vline(xintercept = 0, color = 'black', linetype = 'solid') +
  facet_wrap(~Category, nrow = 1, axes = 'all') +
  labs(x = 'Nascent_log2FC_GB') +
  theme(panel.grid.major = element_line())

ggplot(subset(df_to_plot, Category == 'nsPI_upGB'), aes(x = DEseq_log2FC_GB, y = RNAseq_log2FC)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE, color = 'red', linewidth = 0.5, linetype = 'dashed') +
  geom_hline(yintercept = 1, color = 'blue', linetype = 'dashed') +
  geom_vline(xintercept = 1, color = 'blue', linetype = 'dashed') +
  labs(x = 'Nascent_log2FC_GB') +
  theme(panel.grid.major = element_line())

library(ggvenn)

Nascent_fullGene_DE <- read.xlsx("DiffExpressedGenes_fullGene.xlsx", rowNames = TRUE)
RNAseq_E16_DEs <- read.xlsx('/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/RNAseq_E16_genotypeOnly_noH614_DEgenes_annotated.xlsx', 
                            sheet = 1) %>% distinct(ensembl, .keep_all = TRUE) %>%
  select(c(ensembl, log2FoldChange, padj)) %>% rename(GeneID = ensembl)

RNAseq_UpGenes <- subset(RNAseq_E16_DEs, log2FoldChange > 1 & padj < 0.05 & padj > 1e-50, select = GeneID, drop = TRUE)
Nascent_UpGenes <- subset(Nascent_fullGene_DE, log2FoldChange > 1 & padj < 0.05, select = Ensembl, drop = TRUE)
ggvenn(list('RNA-seq' = RNAseq_UpGenes, 'Nascent' = Nascent_UpGenes), auto_scale = TRUE)


#################### GO Analysis ##########
library(clusterProfiler)
library(org.Mm.eg.db)

wb <- createWorkbook()
addWorksheet(wb, "Nascent_fullGeneUp")
go_DE <- enrichGO(Nascent_UpGenes, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "ENSEMBL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go_DE, showCategory = 10, font.size = 9) + labs(title = 'Nascent_up')
writeData(wb, sheet = 'Nascent_fullGeneUp', as.data.frame(go_DE), rowNames = TRUE)

GO_Plots <- list()
for (name in c('upPI_nsGB', 'nsPI_upGB', 'downPI_upGB', 'downPI_nsGB')){
  cat('Processing:', name, '\n')
  go_DE <- enrichGO(Nascent_GOI[[name]], OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "ENSEMBL", pvalueCutoff = 0.05, readable = TRUE)
  GO_Plots[[name]] <- barplot(go_DE, showCategory = 10, font.size = 9) + labs(title = name)
  addWorksheet(wb, name)
  writeData(wb, sheet = name, as.data.frame(go_DE), rowNames = TRUE)
}

commonDEs <- intersect(RNAseq_UpGenes, Nascent_UpGenes)
go_DE <- enrichGO(commonDEs, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "ENSEMBL", pvalueCutoff = 0.05, readable = TRUE)
addWorksheet(wb, "Shared_withBulkRNA")
writeData(wb, sheet = "Shared_withBulkRNA", as.data.frame(go_DE), rowNames = TRUE)

saveWorkbook(wb, file = 'GO_Nascent_DEgenes.xlsx', overwrite = TRUE)
rm(wb, go_DE, name)
rm(GO_Plots)

#### Plot GO results
go_table <- read.xlsx("GO_Nascent_DEgenes.xlsx", sheet = "Shared_withBulkRNA", rowNames = TRUE)
GO_terms_to_plot <- c("GO:0051321", "GO:0045143", "GO:0140013", "GO:0048477",
                      "GO:0007281", "GO:0007028", "GO:0010526", "GO:0006310")

df_to_plot <- go_table %>% arrange(Count) %>% mutate(Description = factor(Description, levels = Description)) %>%
  mutate(LogPadjust = -log10(p.adjust)) %>%
  subset(ID %in% GO_terms_to_plot)

ggplot(df_to_plot, aes(x = Description, y = Count, fill = LogPadjust)) +
  geom_bar(stat = 'identity', color = 'black') +
  coord_flip() +
  scale_fill_gradient(low = 'white', high = 'darkgreen', limits = c(0, 8)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(x = NULL, y = "Gene count", fill = "-log10(p.adjust)", title = 'Shared Up genes') +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12))
