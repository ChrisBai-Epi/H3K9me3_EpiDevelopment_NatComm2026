library(tidyverse)
library(openxlsx)
library(GenomicFeatures)
library(clusterProfiler)
library(org.Mm.eg.db)

##### Obtain gene lists for C1-C4 methylated clusters and perform GO analysis
Methyl_Clusters <- read.xlsx("DataFiles/Methyl_Annotation/Methyl_clusters.xlsx", colNames = TRUE)
table(Methyl_Clusters$Cluster)
info <- select(Ens_mm39_110, keys = unlist(all_genes$gene_id), keytype = 'GENEID', 
               columns = c('GENEBIOTYPE', 'GENEID', 'GENENAME', 'SYMBOL', 'DESCRIPTION'))
Methyl_Clusters <- left_join(Methyl_Clusters, info, by = join_by(Ensembl == GENEID))
Methyl_Clusters <- dplyr::select(Methyl_Clusters, !c(SYMBOL, GENENAME))

methyl_associated_genes <- list()
for (cluster in c('C1', 'C2', 'C3', 'C4')){
  df <- subset(Methyl_Clusters, Cluster == cluster & GENEBIOTYPE == 'protein_coding')
  methyl_associated_genes[[cluster]] <- df$Symbol
}
rm(df, cluster)

GO_tables <- lapply(methyl_associated_genes, function(glist){
  go <- enrichGO(glist, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
  table <- as.data.frame(go)
})
GO_plots <- lapply(methyl_associated_genes, function(glist){
  go <- enrichGO(glist, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
  barplot(go)
})

write.xlsx(GO_tables[['C1']], "DataFiles/GOAnalyses/GO_methyl_C1_genes.xlsx", rowNames = TRUE)
write.xlsx(GO_tables[['C2']], "DataFiles/GOAnalyses/GO_methyl_C2_genes.xlsx", rowNames = TRUE)
write.xlsx(GO_tables[['C3']], "DataFiles/GOAnalyses/GO_methyl_C3_genes.xlsx", rowNames = TRUE)
write.xlsx(GO_tables[['C4']], "DataFiles/GOAnalyses/GO_methyl_C4_genes.xlsx", rowNames = TRUE)

## Combine C1~C4 genes and perform GO anlaysis
combinedList <- Reduce(union, methyl_associated_genes)
go <- enrichGO(combinedList, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go, showCategory = 15)

## Check the expression of these methylation-associated genes
RNAseq_E16_dataAll <- read.xlsx('/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/RNAseq_E16_genotypeOnly_noH614_DEgenes_annotated.xlsx', 
                            sheet = 1)
df_to_plot <- RNAseq_E16_dataAll %>% subset(grepl('Zfp', GeneID)) %>%
  mutate(Label = if_else(log2FoldChange > 1, GeneID, NA))
df_to_plot <- RNAseq_E16_dataAll %>% subset(GeneID %in% methyl_associated_genes$C1) %>%
  mutate(Label = if_else(log2FoldChange > 1, GeneID, NA))

ggplot(df_to_plot, aes(x = log10(baseMean), y = log2FoldChange, label = Label)) +
  geom_point() +
  geom_hline(yintercept = 1, color = 'blue', linetype = 'dashed') +
  geom_hline(yintercept = -1, color = 'blue', linetype = 'dashed') +
  geom_text(nudge_x = 0.2, nudge_y = 0.1, size = 3, color = 'black')


####### Try defining H3K9me3-associated genes as genes that overlap epic2 peaks (including 10kb flanking regions)
## Did not use eventually
all_genes_expanded10kb <- resize(all_genes, fix = 'center', width = width(all_genes) + 20000)
highScore_peaks <- subset(Epic2_peaks, Score > 150)
idx <- findOverlaps(all_genes_expanded10kb, highScore_peaks, select = 'first')
methyl_associated_genes <- subset(all_genes, subset = !is.na(idx) & Biotype == 'protein_coding')

go <- enrichGO(methyl_associated_genes$gene_id, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "ENSEMBL", pvalueCutoff = 0.05, readable = TRUE)
view(as.data.frame(go))


###### Visualization of GO analysis result ######
## Read in saved GO result tables
GO_tables <- list()
for (name in c('C1', 'C2', 'C3', 'C4')){
  path <- paste0("DataFiles/GOAnalyses/GO_methyl_", name, "_genes.xlsx")
  GO_tables[[name]] <- read.xlsx(path)
}
rm(name, path)

## SimplifyEnrichment analysis
library(simplifyEnrichment)

local({
  ids <- subset(as.data.frame(GO_tables$C4), ONTOLOGY == 'BP', select = ID, drop = TRUE)
  mat = GO_similarity(ids)
  simplifyGO(mat)
})

## Select GO terms for dot plot
selectedGO <- c('GO:0001227', 'GO:0120261', 'GO:0019882', 'GO:0004497',
                'GO:0048515', 'GO:0007281', 'GO:0048663', 'GO:0001676', 'GO:0009952',
                'GO:0045088', 'GO:0006068', 'GO:0007389', 'GO:0035107')
selectedGO <- c('GO:0001217', 'GO:0120261', 'GO:0019882', 'GO:0004497', 'GO:0031731')

GO_list_all <- bind_rows(GO_tables[1:4], .id = 'Cluster') %>% subset(ID %in% selectedGO) %>%
  mutate(ID = factor(ID, levels = selectedGO)) %>% arrange(desc(ID)) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

ggplot(GO_list_all, aes(x = Cluster, y = Description, color = qvalue)) +
  geom_point(aes(size = Count)) +
  scale_size(range = c(3,9), breaks = c(4, 7, 10)) +
  scale_color_gradient(low = 'darkgreen', high = 'gray90', limits = c(0, 0.05)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

rm(GO_list_all)

## Select GO terms for barplot
## Barplot does not contain information of which cluster each term is enriched in, so not used.
selectedGO <- c('GO:0001217', 'GO:0120261', 'GO:0019882', 'GO:0004497', 'GO:0031731')

go_table <- Reduce(rbind, GO_tables[1:3])
df_to_plot <- go_table %>% arrange(Count) %>% mutate(Description = factor(Description, levels = Description)) %>%
  mutate(LogPadjust = -log10(p.adjust)) %>%
  subset(ID %in% selectedGO)

ggplot(df_to_plot, aes(x = Description, y = Count, fill = LogPadjust)) +
  geom_bar(stat = 'identity', color = 'black') +
  coord_flip() +
  scale_fill_gradient(low = 'white', high = 'darkgreen', limits = c(0, 5)) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 13)) +
  labs(x = NULL, y = "Gene count", fill = "-log10(p.adjust)") +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12))
