library(DESeq2)
library(org.Mm.eg.db)
library(tidyverse)
library(scuttle)
library(openxlsx)


###### DESeq2 Analysis #####
#Read in data
rawCounts <- as.matrix(read.table("6029D_rawCounts.txt"))
metaInfo <- read.table("6029D_metaData.txt", header = TRUE)

metaInfo$Genotype <- factor(metaInfo$Genotype, levels = c("Ctrl", "TKO"))
metaInfo$Sex <- factor(metaInfo$Sex, levels = c("F", "M"))
metaInfo$Batch <- factor(metaInfo$Batch)

dds_without_H614 <- DESeqDataSetFromMatrix(countData = rawCounts, colData = metaInfo, design = ~ Genotype)
dds_without_H614 <- dds_without_H614[,c(-1,-5)]

#Getting Annotation
anno_symbol <- mapIds(org.Mm.eg.db, keys = row.names(dds_without_H614), keytype = 'ENSEMBL', column = 'SYMBOL')
anno_symbol <- uniquifyFeatureNames(names(anno_symbol), anno_symbol)
rowData(dds_without_H614)$ensembl <- row.names(dds_without_H614)
rowData(dds_without_H614)$symbol <- anno_symbol
row.names(dds_without_H614) <- rowData(dds_without_H614)$symbol

##Analysis
keep <- rowSums(counts(dds_without_H614)) >= 100
dds_without_H614 <- dds_without_H614[keep,]
rm(keep)

dds_without_H614 <- DESeq(dds_without_H614)
res <- results(dds_without_H614, contrast = c('Genotype', 'TKO', 'Ctrl'), alpha = 0.05, saveCols = "ensembl")
#res <- results(dds_without_H614, contrast = c('Sex', 'M', 'F'), alpha = 0.05, saveCols = "ensembl")
resLFC <- lfcShrink(dds_without_H614, coef = "Genotype_TKO_vs_Ctrl", type = 'apeglm')

res$padj <- ifelse(is.na(res$padj), 1, res$padj)
DEgenes <- res[res$padj<0.05 & abs(res$log2FoldChange) > 1, ]

write_DEseq_result(res, 'E16_genotypeOnly_noH614_DEgenes.xlsx')
addGeneInfo('E16_genotypeOnly_noH614_DEgenes.xlsx', ensembl = Ens_mm39_110, ID_type = 'SYMBOL')
file.remove('E16_genotypeOnly_noH614_DEgenes.xlsx')


##### Volcano Plot #####
library(EnhancedVolcano)

Gene_of_interest = c('Col28a1', 'Kyat3', 'Dpep1', 'Tcfl5', 'Ly6e')

res_to_plot <- subset(as.data.frame(res), padj > 1e-50)
down_genes_color = ifelse(res_to_plot$padj > 0.05 | abs(res_to_plot$log2FoldChange) <= 1, regColor[['ns']], ifelse(res_to_plot$log2FoldChange > 1, regColor[['up']], regColor[['down']]))
down_genes_color[is.na(down_genes_color)] <- regColor[['ns']]
names(down_genes_color)[down_genes_color == regColor[['up']]] <- 'Up regulated'
names(down_genes_color)[down_genes_color == regColor[['down']]] <- 'Down regulated'
names(down_genes_color)[down_genes_color == regColor[['ns']]] <- 'Not significant'

EnhancedVolcano(res_to_plot, lab = row.names(res_to_plot), x = 'log2FoldChange', y = 'pvalue', title = 'TKO vs Ctrl', FCcutoff = 1, pCutoff = 0.05, pCutoffCol = 'padj', 
                colCustom = down_genes_color, legendPosition = 'none', selectLab = Gene_of_interest,
                labSize = 5, pointSize = 3, axisLabSize = 20, subtitle = NULL, caption = bquote(~Log[2]~ "fold change cutoff, 1; adjusted p-value cutoff, 0.05"))

#Draw PCA plots
vsd <- vst(dds_without_H614, blind = FALSE)
pca <- plotPCA(vsd, pcsToUse = c(1,2), intgroup = 'Genotype', ntop = 2000, returnData = FALSE)
pca +
  scale_color_manual(values = genoColor) +
  labs(colour = 'Genotype') +
  theme(panel.grid.major = element_line())

save.image()

###### GO analysis. ######
library(clusterProfiler)

RNAseq_E16_DEs <- read.xlsx('/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/RNAseq_E16_genotypeOnly_noH614_DEgenes_annotated.xlsx', 
                            sheet = 2)
RNAseq_E16_upGenes <- subset(RNAseq_E16_DEs, log2FoldChange > 1, select = 'GeneID', drop = TRUE)

go_up <- enrichGO(RNAseq_E16_upGenes, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
go_table <- as.data.frame(go_up)
write.xlsx(go_table, file = "GO_UpGenes.xlsx")

#GO_terms_to_plot <- c("GO:0051321", "GO:0060047", "GO:0006936", "GO:0098926", "GO:0030239", "GO:0007281", "GO:1901605", "GO:0032197", "GO:0098984", "GO:1901605")
#As commented by Doina, remove retrotransposition term
GO_terms_to_plot <- c("GO:0051321", "GO:0060047", "GO:0006936", "GO:0098926", "GO:0030239", "GO:0007281", "GO:1901605", "GO:0098984", "GO:1901605")
df_to_plot <- go_table %>% arrange(Count) %>% mutate(Description = factor(Description, levels = Description)) %>%
  mutate(LogPadjust = -log10(p.adjust)) %>%
  subset(ID %in% GO_terms_to_plot)
ggplot(df_to_plot, aes(x = Description, y = Count, fill = LogPadjust)) +
  geom_bar(stat = 'identity', color = 'black') +
  coord_flip() +
  scale_fill_gradient(low = 'white', high = 'darkgreen', limits = c(0, 5)) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 13)) +
  labs(x = NULL, y = "Gene count", fill = "-log10(p.adjust)") +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12))
