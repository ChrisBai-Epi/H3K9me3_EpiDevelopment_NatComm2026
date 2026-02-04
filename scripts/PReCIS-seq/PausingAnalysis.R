library(tidyverse)
library(DESeq2)
library(openxlsx)

######PI FC analysis using DESeq2
all(row.names(GeneBodyCounts) == row.names(PauseCounts))

Counts <- cbind(GeneBodyCounts[,c(3,4,7,8)], PauseCounts[,c(3,4,7,8)])
colnames(Counts) <- c('GB_CT_IP_A', 'GB_CT_IP_B', 'GB_KO_IP_A', 'GB_KO_IP_B',
                      'PP_CT_IP_A', 'PP_CT_IP_B', 'PP_KO_IP_A', 'PP_KO_IP_B')

metaInfo <- data.frame(Name = names(Counts),
                       Genotype = rep(c('Ctrl', 'Ctrl', 'TKO', 'TKO'), 2),
                       Region = rep(c('GB', 'PR'), each = 4))
metaInfo$Sample <- paste0(metaInfo$Genotype, '_', metaInfo$Region)
metaInfo$Sample <- factor(metaInfo$Sample, levels = unique(metaInfo$Sample))

dds_for_PI <- DESeqDataSetFromMatrix(countData = Counts, colData = metaInfo, design = ~ Sample)
rowData(dds_for_PI)$Ensembl <- names(dds_for_PI)
info <- AnnotationDbi::select(Ens_mm39_110, keys = rowData(dds_for_PI)$Ensembl, keytype = 'GENEID', columns = c('GENEID', 'SYMBOL', 'GENEBIOTYPE', 'DESCRIPTION'))
rowData(dds_for_PI)$Symbol <- info$SYMBOL

#Define NormalizationFactor, not SizeFactors (not working since PR and GB region have different length)
size_factor <- readCountSummary$exp_reads[c(3,4,7,8)]
size_factor <- size_factor/exp(mean(log(size_factor)))
GB_norm_factor <- outer(gene_length, size_factor)
PR_norm_factor <- outer(rep(0.25, length(gene_length)), size_factor)
normFactors <- cbind(GB_norm_factor, PR_norm_factor)
normFactors <- normFactors / exp(rowMeans(log(normFactors)))
normalizationFactors(dds_for_PI) <- normFactors
rm(normFactors)

dds_for_PI <- DESeq(dds_for_PI)
dds <- dds_for_PI[rowSums(counts(dds_for_PI))>=80,]
#Formula should be KO_PI/CT_PI, which is KO_PR*CT_GB/(KO_GB*CT_PR), but intercept is base level and all other coefficents
#are changes, so log2(CT_PR) - log2(CT_GB) is actually Beta_CT_PR, which is "Sample_Ctrl_PR_vs_Ctrl_GB, so the contrast
#should be c(0, -1, -1, 1), while 0 represent intercept shouldn't be included for calculating changes
res <- results(dds, contrast = c(0, -1, -1, 1), alpha = 0.05, saveCols = c("Ensembl", 'Symbol'))
write_DEseq_result(res, filename = 'Diff_PI_genes_byDESeq2.xlsx')

res <- as.data.frame(res)
res$padj <- ifelse(is.na(res$padj), 1, res$padj)
UpPI_genes <- subset(res, pvalue < 0.01 & log2FoldChange > 1)
DownPI_genes <- subset(res, pvalue < 0.01 & log2FoldChange < -1)


###Scatter plots for PI vs GB
Merged_Count_Table <- data.frame(row.names = rownames(PauseCounts),
                                 GB_CT_IP = GeneBodyCounts$CT_IP_A + GeneBodyCounts$CT_IP_B,
                                 GB_KO_IP = GeneBodyCounts$KO_IP_A + GeneBodyCounts$KO_IP_B,
                                 PS_CT_IP = PauseCounts$CT_IP_A + PauseCounts$CT_IP_B,
                                 PS_KO_IP = PauseCounts$KO_IP_A + PauseCounts$KO_IP_B)

###Calculate RPKM and then PI and then log2FC
Norm_factor <- c(0.03232062, 0.09453583) #This is calculated in a seperate excel file for spike-in normalization
Merged_PI_Table <- Merged_Count_Table %>% mutate(GB_CT_IP = GB_CT_IP/gene_length*Norm_factor[1], 
                                                 GB_KO_IP = GB_KO_IP/gene_length*Norm_factor[2],
                                                 PS_CT_IP = PS_CT_IP*1000/250*Norm_factor[1], 
                                                 PS_KO_IP = PS_KO_IP*1000/250*Norm_factor[2]) %>%
  filter(rowSums(Merged_Count_Table) >= 80) %>% mutate(PI_CT = (PS_CT_IP + 0.01)/(GB_CT_IP + 0.01),
                                           PI_KO = (PS_KO_IP + 0.01)/(GB_KO_IP + 0.01)) %>%
  mutate(PI_log2FC = log2(PI_KO/PI_CT),
         GB_log2FC = log2((GB_KO_IP + 0.01)/(GB_CT_IP + 0.01)))

ggplot(Merged_PI_Table, aes(x = GB_log2FC, y = PI_log2FC)) +
  geom_point(size = 1)+
  coord_fixed(ratio = 1)+
  geom_vline(xintercept = 1, linetype = 'dotted', color = 'red') +
  geom_vline(xintercept = -1, linetype = 'dotted', color = 'red') +
  geom_hline(yintercept = 1, linetype = 'dotted', color = 'red') +
  geom_hline(yintercept = -1, linetype = 'dotted', color = 'red') +
  labs(x="Log2(GB-KO/GB-CT)",
       y="Log2(PI-KO/PI-CT)",
       title = "By Spike-in")

### Obtain a list of unexpressed genes
Unexpressed <- rownames(Merged_Count_Table[rowSums(Merged_Count_Table)<=10,])
saveRDS(Unexpressed, file = 'Unexpressed_GeneList.rds')
rm(Unexpressed)

#####Scatter plot2, using statiscally tested log2FC data
PI_Res_DESeq <- read.xlsx("Diff_PI_genes_byDESeq2.xlsx", sheet = 1, rowNames = TRUE) %>%
  mutate(padj = ifelse(is.na(padj), 1, padj))
DE_Res_SpikeIn <- read.xlsx("DiffExpressedGenes_SpikeInNorm_annotated.xlsx", sheet = 1, rowNames = TRUE) %>%
  mutate(padj = ifelse(is.na(padj), 1, padj))

highPPCGenes <- rownames(subset(Merged_Count_Table, PS_CT_IP > 5 | PS_KO_IP > 5)) # No need to filter

FC_table_w_test <- merge(Merged_PI_Table[,7:8], PI_Res_DESeq[,c(2,5)], by = 'row.names', sort = FALSE) %>%
  column_to_rownames("Row.names") %>% merge(DE_Res_SpikeIn[,c(2,6)], by = 'row.names', sort = FALSE) %>%
  rename(DEseq_log2FC_PI = log2FoldChange.x, DEseq_log2FC_GB = log2FoldChange.y) %>%
  column_to_rownames("Row.names") %>%
  mutate(Reg_PI = if_else(pvalue > 0.01 | abs(DEseq_log2FC_PI) < 1, 'ns', if_else(DEseq_log2FC_PI > 1, 'up', 'down')),
         Reg_GB = if_else(padj > 0.05 | abs(DEseq_log2FC_GB) < 1, 'ns', if_else(DEseq_log2FC_GB > 1, 'up', 'down')))

FC_table_w_test <- FC_table_w_test %>% 
  mutate(Reg_GB = factor(Reg_GB, levels = c('ns', 'down', 'up'))) %>%
  mutate(Reg_PI = factor(Reg_PI, levels = c('ns', 'down', 'up'))) %>%
  arrange(Reg_GB, Reg_PI)

ggplot(FC_table_w_test, aes(x = DEseq_log2FC_GB, y = DEseq_log2FC_PI, colour = Reg_PI, fill = Reg_GB)) +
  geom_point(size = 2, shape = 21) +
  scale_fill_manual(values = c('gray', 'lightblue', 'orange')) +
  scale_color_manual(values = c('gray', 'black', 'red')) +
  coord_fixed(ratio = 1) +
  geom_vline(xintercept = 1, linetype = 'dotted', color = 'black') +
  geom_vline(xintercept = -1, linetype = 'dotted', color = 'black') +
  geom_hline(yintercept = 1, linetype = 'dotted', color = 'black') +
  geom_hline(yintercept = -1, linetype = 'dotted', color = 'black') +
  labs(x = 'Log2(TKO/Ctrl) for GB',
       y = 'Log2(TKO/Ctrl) for PI',
       fill = 'GB', color = 'PI')


info <- column_to_rownames(info, "GENEID")
FC_table_w_test <- add_column(FC_table_w_test, Symbol = info[rownames(FC_table_w_test), 1], .before = 1)
FC_table_w_test <- FC_table_w_test %>%
  mutate(Category = paste0(Reg_PI, 'PI_', Reg_GB, 'GB')) %>%
  mutate(Category = if_else(Category %in% c('nsPI_nsGB', 'nsPI_downGB'), 'others', Category))


write.xlsx(FC_table_w_test, file = "FoldChangeTable.xlsx", rowNames = TRUE)
addGeneInfo(filename = "FoldChangeTable.xlsx", Ens_mm39_110, ID_type = 'GENEID')

