library(DESeq2)
library(openxlsx)
library(tidyverse)


#Get counts from the full gene region with Importing_and_Processing.R

metaInfo <- data.frame(Sample = names(GeneBodyCounts),
                       Genotype = c(rep('Ctrl', 4), rep('TKO', 4)),
                       Type = rep(c('Input', 'Input', 'IP', 'IP'), 2))
metaInfo$Genotype <- factor(metaInfo$Genotype, levels = c('Ctrl', 'TKO'))
metaInfo$Type <- factor(metaInfo$Type, levels = c('Input', 'IP'))

counts <- FullCounts[, c(3,4,7,8)]
dds_full_counts <- DESeqDataSetFromMatrix(countData = counts, colData = metaInfo[c(3,4,7,8),], design = ~ Genotype)
rowData(dds_full_counts)$Ensembl <- names(dds_full_counts)
rm(counts)

## Normalize with Spike-in
dds_full_spikeIn_Norm <- dds_full_counts
normFactors <- outer(gene_full_length, 1/NF_SpikeIn[c(3,4,7,8)])
normFactors <- normFactors / exp(rowMeans(log(normFactors)))
normalizationFactors(dds_full_spikeIn_Norm) <- normFactors
dds_full_spikeIn_Norm <- DESeq(dds_full_spikeIn_Norm)
rm(normFactors)

dds <- dds_full_spikeIn_Norm[rowSums(counts(dds_full_spikeIn_Norm))>=80,]
rowData(dds) <- cbind(rowData(dds), info[rownames(dds),])
res <- results(dds, contrast = c('Genotype', 'TKO', 'Ctrl'), alpha = 0.05,
               saveCols = c("Ensembl", "SYMBOL", "GENEBIOTYPE", "DESCRIPTION"))

write.xlsx(as.data.frame(res), file = "DiffExpressedGenes_fullGene.xlsx", rowNames = TRUE)

## Plot Volcano plot
library(EnhancedVolcano)

Nascent_fullGene_DE <- read.xlsx("DiffExpressedGenes_fullGene.xlsx", rowNames = TRUE)
Gene_of_interest = c('Col28a1', 'Kyat3', 'Dpep1', 'Tcfl5', 'Stag3', 'Ly6e', 'Nat8f2')

df_to_plot <- as.data.frame(Nascent_fullGene_DE)
regColor <- c('up' = "red", 'ns' = "grey", 'down' = "blue")
genes_color = ifelse(df_to_plot$padj > 0.05 | abs(df_to_plot$log2FoldChange) <= 1, regColor[['ns']], ifelse(df_to_plot$log2FoldChange > 1, regColor[['up']], regColor[['down']]))
genes_color[is.na(genes_color)] <- regColor[['ns']]
names(genes_color)[genes_color == regColor[['up']]] <- 'Up regulated'
names(genes_color)[genes_color == regColor[['down']]] <- 'Down regulated'
names(genes_color)[genes_color == regColor[['ns']]] <- 'Not significant'

EnhancedVolcano(df_to_plot, lab = df_to_plot$SYMBOL, x = 'log2FoldChange', y = 'pvalue', title = 'TKO vs Ctrl', FCcutoff = 1, pCutoff = 0.05, pCutoffCol = 'padj', 
                colCustom = genes_color, legendPosition = 'none', selectLab = Gene_of_interest,
                labSize = 5, pointSize = 3, axisLabSize = 20, subtitle = NULL, caption = NULL,
                raster = TRUE)







