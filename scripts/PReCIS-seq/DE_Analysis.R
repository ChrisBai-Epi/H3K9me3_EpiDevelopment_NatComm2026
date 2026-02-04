library(org.Mm.eg.db)
library(DESeq2)

FoldChange_table <- data.frame(ID = Genes_Granges$gene_id) %>% column_to_rownames('ID')
FoldChange_table$Col_1 <- NA
FoldChange_table$Col_2 <- NA
#These are just placeholders so that the following codes can run smoothly

metaInfo <- data.frame(Sample = names(GeneBodyCounts),
                       Genotype = c(rep('Ctrl', 4), rep('TKO', 4)),
                       Type = rep(c('Input', 'Input', 'IP', 'IP'), 2))
metaInfo$Genotype <- factor(metaInfo$Genotype, levels = c('Ctrl', 'TKO'))
metaInfo$Type <- factor(metaInfo$Type, levels = c('Input', 'IP'))


####Using BRgenomics count table and DESeq2
counts <- GeneBodyCounts[, c(3,4,7,8)]
dds_GB_counts <- DESeqDataSetFromMatrix(countData = counts, colData = metaInfo[c(3,4,7,8),], design = ~ Genotype)
rowData(dds_GB_counts)$Ensembl <- names(dds_GB_counts)
rm(counts)

#Auto-normalization with DEseq default method
dds_GB_autoNorm <- DESeq(dds_GB_counts)
dds <- dds_GB_autoNorm[rowSums(counts(dds_GB_autoNorm))>=80,]
res <- results(dds, contrast = c('Genotype', 'TKO', 'Ctrl'), alpha = 0.05, saveCols = "Ensembl")
FoldChange_table[res$Ensembl, 1] <- res$log2FoldChange
colnames(FoldChange_table)[1] <- 'BRG_DESeq_RPM'

#Use calculated NF_RPM for RPM normalization
dds_GB_RPM_Norm <- dds_GB_counts
normFactors <- outer(gene_length, 1/NF_RPM[c(3,4,7,8)])
normFactors <- normFactors / exp(rowMeans(log(normFactors)))
normalizationFactors(dds_GB_RPM_Norm) <- normFactors
dds_GB_RPM_Norm <- DESeq(dds_GB_RPM_Norm)
rm(normFactors)

dds <- dds_GB_RPM_Norm[rowSums(counts(dds_GB_RPM_Norm))>=80,]
res <- results(dds, contrast = c('Genotype', 'TKO', 'Ctrl'), alpha = 0.05, saveCols = "Ensembl")
FoldChange_table[res$Ensembl, 1] <- res$log2FoldChange
colnames(FoldChange_table)[1] <- 'BRG_DESeq_RPM'

###Using spikeIn normalizing matrix
dds_GB_spikeIn_Norm <- dds_GB_counts
normFactors <- outer(gene_length, 1/NF_SpikeIn[c(3,4,7,8)])
normFactors <- normFactors / exp(rowMeans(log(normFactors)))
normalizationFactors(dds_GB_spikeIn_Norm) <- normFactors
dds_GB_spikeIn_Norm <- DESeq(dds_GB_spikeIn_Norm)
rm(normFactors)

dds <- dds_GB_spikeIn_Norm[rowSums(counts(dds_GB_spikeIn_Norm))>=80,]
res <- results(dds, contrast = c('Genotype', 'TKO', 'Ctrl'), alpha = 0.05, saveCols = "Ensembl")
FoldChange_table[res$Ensembl, 2] <- res$log2FoldChange
colnames(FoldChange_table)[2] <- 'BRG_DESeq_SpikeIn'


####Using DESeq2###################
view(StarCounts_Mouse)
view(StarCounts_Fly)

counts <- StarCounts_Mouse[,c(4, 5, 8, 9)]
rownames(counts) <- StarCounts_Mouse$ID
dds_mouse <- DESeqDataSetFromMatrix(countData = counts, colData = metaInfo[c(3,4,7,8),], design = ~ Genotype)
rm(counts)

rowData(dds_mouse)$Ensembl <- names(dds_mouse)

dds_mouse <- DESeq(dds_mouse)

dds <- dds_mouse[rowSums(counts(dds_mouse))>=80,]
res <- results(dds, contrast = c('Genotype', 'TKO', 'Ctrl'), alpha = 0.05, saveCols = "Ensembl")
FoldChange_table[res$Ensembl, 3] <- res$log2FoldChange
colnames(FoldChange_table)[3] <- 'Star_DESeq_RPM'

###Using DESeq2 SpikeIn
counts <- StarCounts[,c(4, 5, 8, 9)]
rownames(counts) <- StarCounts$ID
dds_star_all <- DESeqDataSetFromMatrix(countData = counts, colData = metaInfo[c(3,4,7,8),], design = ~ Genotype)
rowData(dds_star_all)$Ensembl <- names(dds_star_all)
rm(counts)

dds_star_all <- estimateSizeFactors(dds_star_all, controlGenes = grepl("^FBgn", names(dds_star_all)))
dds_star_all <- DESeq(dds_star_all)

dds <- dds_star_all[rowSums(counts(dds_star_all))>=80,]
res <- results(dds, contrast = c('Genotype', 'TKO', 'Ctrl'), alpha = 0.05, saveCols = "Ensembl")
FoldChange_table[res$Ensembl, 4] <- res$log2FoldChange
colnames(FoldChange_table)[4] <- 'Star_DESeq_SpikeIn'
FoldChange_table <- FoldChange_table[grep("^ENSMUSG", rownames(FoldChange_table)),]


#Direct fold change calculation from count matrix (Not using anymore)
GB_FC <- GeneBodyCounts_RPKM %>% mutate(FC = (KO_IP_A + KO_IP_B + 0.1) / (CT_IP_A + CT_IP_B + 0.1)) %>% mutate(Log2FC = log2(FC))
GB_FC$Log2FC[rowSums(GeneBodyCounts) <= 80] <- NA
FoldChange_table[rownames(GB_FC), 5] <- GB_FC$Log2FC

GB_FC <- GeneBodyCounts_SPKM %>% mutate(FC = (KO_IP_A + KO_IP_B + 0.1) / (CT_IP_A + CT_IP_B + 0.1)) %>% mutate(Log2FC = log2(FC))
GB_FC$Log2FC[rowSums(GeneBodyCounts) <= 80] <- NA
FoldChange_table[rownames(GB_FC), 6] <- GB_FC$Log2FC

colnames(FoldChange_table)[5:6] <- c('GB_RPKM', 'GB_SPKK')

#Direct fold change calculation from merged count matrix
Merged_GB_IP_Counts <- data.frame(ID = Genes_Granges$gene_id) %>% column_to_rownames('ID')

Merged_GB_IP_Counts$CT_IP_Counts <- GeneBodyCounts$CT_IP_A + GeneBodyCounts$CT_IP_B
Merged_GB_IP_Counts$KO_IP_Counts <- GeneBodyCounts$KO_IP_A + GeneBodyCounts$KO_IP_B

ReadCounts_for_merged <- column_to_rownames(readCountSummary, "sample") %>% t() %>% as.data.frame() %>%
  mutate(CT_IP_merged = CT_IP_A + CT_IP_B, KO_IP_merged = KO_IP_A + KO_IP_B)

NF_for_merged <- as.matrix(c(1e6, 1e3) / ReadCounts_for_merged[2:3, 9:10])

normalizedCounts_RPKM <- Merged_GB_IP_Counts * outer(1/gene_length, NF_for_merged['exp_reads',])
normalizedCounts_SPKK <- Merged_GB_IP_Counts * outer(1/gene_length, NF_for_merged['spike_reads',])

colnames(normalizedCounts_RPKM) <- c('CT_IP_RPKM', 'KO_IP_RPKM')
colnames(normalizedCounts_SPKK) <- c('CT_IP_SPKK', 'KO_IP_SPKK')

Merged_GB_IP_allData <- Merged_GB_IP_Counts %>% cbind(normalizedCounts_RPKM) %>% cbind(normalizedCounts_SPKK) %>%
  mutate(FC_RPKM = log2((KO_IP_RPKM + 0.01) / (CT_IP_RPKM + 0.01)), FC_SPKK = log2((KO_IP_SPKK + 0.01) / (CT_IP_SPKK + 0.01)))
Merged_GB_IP_allData$FC_RPKM[rowSums(Merged_GB_IP_Counts) <= 80] <- NA
Merged_GB_IP_allData$FC_SPKK[rowSums(Merged_GB_IP_Counts) <= 80] <- NA

FoldChange_table[rownames(Merged_GB_IP_allData), 5] <- Merged_GB_IP_allData$FC_RPKM
FoldChange_table[rownames(Merged_GB_IP_allData), 6] <- Merged_GB_IP_allData$FC_SPKK
colnames(FoldChange_table)[5:6] <- c('GB_RPKM', 'GB_SPKK')


######Compare Different Methods#################

#Remove all NA rows
FC_table <- FoldChange_table %>% filter(!if_any(everything(), is.na))

#Correlation matrix
cor_mat <- cor(FC_table, method = "pearson")
ggcorrplot::ggcorrplot(cor_mat, method = "circle")

#Scatter plots
library(ggpmisc)
library(rlang)

ggplot(FC_table, aes(x = BRG_DESeq_RPM, y = BRG_DESeq_SpikeIn)) +
  geom_point() +
  coord_fixed(ratio = 1) +
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 'dashed') +
  geom_smooth(method = 'lm', se = FALSE, color = 'blue', linewidth = 0.5) +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = '~~~')), formula = y ~ x, parse = TRUE, color = 'black') +
  labs(x = "Log2(TKO/Ctrl) by total reads (RPKM)", y = "Log2(TKO/Ctrl) by Spike-in (SPKK)", title = "DESeq2") +
  theme(panel.grid.major = element_line(),
        panel.grid.minor = element_line())

plot_correlation_CB <- function(df, x_val, y_val){
  require(ggpmisc)
  ggplot(FC_table, aes(x = {{ x_val }}, y = {{ y_val }})) +
    geom_point() +
    coord_fixed(ratio = 1)+
    geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 'dashed') +
    geom_smooth(method = 'lm', se = FALSE, color = 'blue', linewidth = 0.5) +
    stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = '~~~')), formula = y ~ x, parse = TRUE, color = 'black') +
    theme(panel.grid.major = element_line(),
          panel.grid.minor = element_line())
}

plot_correlation_CB(FC_table, BRG_DESeq_RPM, BRG_DESeq_SpikeIn) +
  labs(x = "Log2(TKO/Ctrl) by total reads (RPKM)", y = "Log2(TKO/Ctrl) by Spike-in (SPKK)", title = "DESeq2")
plot_correlation_CB(FC_table, GB_RPKM, GB_SPKK) +
  labs(x = "Log2(TKO/Ctrl) by RPKM", y = "Log2(TKO/Ctrl) by Spike-in norm", title = "Manual")
plot_correlation_CB(FC_table, BRG_DESeq_RPM, GB_RPKM)
plot_correlation_CB(FC_table, BRG_DESeq_SpikeIn, GB_SPKK)

#####Check DE genes
dds <- dds_GB_spikeIn_Norm[rowSums(counts(dds_GB_spikeIn_Norm))>=80,]
DE_file_name <- "DiffExpressedGenes_SpikeInNorm.xlsx"

dds <- dds_GB_RPM_Norm[rowSums(counts(dds_GB_RPM_Norm))>=80,]
DE_file_name <- "DiffExpressedGenes_RPKM.xlsx"

info <- AnnotationDbi::select(Ens_mm39_110, keys = rowData(dds)$Ensembl, keytype = 'GENEID', columns = c('GENEID', 'SYMBOL', 'GENEBIOTYPE', 'DESCRIPTION'))
rowData(dds)$Symbol <- info$SYMBOL
res <- results(dds, contrast = c('Genotype', 'TKO', 'Ctrl'), alpha = 0.05, saveCols = c("Ensembl", "Symbol"))

write_DEseq_result(res, DE_file_name)
addGeneInfo(DE_file_name, Ens_mm39_110, 'GENEID')

rm(DE_file_name)

