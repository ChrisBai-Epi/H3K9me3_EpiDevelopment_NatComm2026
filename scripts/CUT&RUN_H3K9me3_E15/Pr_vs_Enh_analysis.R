library(tidyverse)
library(openxlsx)
library(GenomicRanges)
library(clusterProfiler)
library(org.Mm.eg.db)

### This script aims to explore genes directly repressed by H3K9me3 at promoters vs those by enhancers

##### With RNA-seq data #####
RNAseq_E16_DEs <- read.xlsx('/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/RNAseq_E16_genotypeOnly_noH614_DEgenes_annotated.xlsx',
                            sheet = 1)
Peak_distance <- read.table("DataFiles/Methyl_Annotation/ClosestEpic2Peak_toPromoters.bed", sep = "\t") %>%
  select(c(V4, V10, V11, V13)) %>% rename(ensembl = V4, PeakScore = V11, Distance = V13, PeakID = V10) %>% filter(PeakScore > 0) %>%
  merge(RNAseq_E16_DEs, ., by = 'ensembl', all = FALSE) %>%
  mutate(Group = if_else(log2FoldChange > 1 & padj < 0.05, "Upregulated", "Others")) %>%
  mutate(Group = factor(Group, levels = c('Upregulated', 'Others'))) %>%
  mutate(LogDistance = case_when(Distance > 0 ~ log10(Distance),
                                 Distance == 0 ~ 0,
                                 Distance < 0 ~ -log10(-Distance))) %>%
  mutate(Reg_element = if_else(abs(Distance) <= 500, 'Promoter', 'Enhancer')) %>%
  arrange(desc(Group))

Upregulated_genes <- Peak_distance %>% subset(Group == 'Upregulated')

info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'description'), mart = enbl, filters = 'ensembl_gene_id', values = Upregulated_genes$ensembl, uniqueRows = TRUE)
Upregulated_genes <- subset(Upregulated_genes, !duplicated(ensembl))
Upregulated_genes$GeneID <- info$external_gene_name
Upregulated_genes$GENEBIOTYPE <- info$gene_biotype
Upregulated_genes$DESCRIPTION <- info$description
Upregulated_genes$GENEID <- Upregulated_genes$ensembl
write.xlsx(Upregulated_genes, file = 'DataFiles/GeneAnalysis_Promoter_Enhancer/BulkRNAseq_Upreg_genes.xlsx')

#GO analysis yeilds no terms for promoter regulated genes yields no terms
Upregulated_genes <- read.xlsx("DataFiles/GeneAnalysis_Promoter_Enhancer/BulkRNAseq_Upreg_genes.xlsx")
promoter_reg_genes <- subset(Upregulated_genes, Reg_element == 'Promoter' & GENEBIOTYPE == 'protein_coding')

enhancer_reg_genes <- subset(Upregulated_genes, Reg_element == 'Enhancer')
go <- enrichGO(enhancer_reg_genes$ensembl, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "ENSEMBL", pvalueCutoff = 0.05, readable = TRUE)
barplot(go, showCategory = 15)
write.xlsx(as.data.frame(go), file = 'DataFiles/GeneAnalysis_Promoter_Enhancer/GO_BulkRNAseq_Enhancer_reg_genes.xlsx')

rm(Upregulated_genes, enhancer_reg_genes, promoter_reg_genes)

#Refine enhancer-regulation with C4 cluster signature
Methyl_Cluster_info <- read.xlsx("DataFiles/Methyl_Annotation/Methyl_clusters.xlsx", colNames = TRUE) %>%
  select(Ensembl, Cluster)
enhancer_reg_genes <- Upregulated_genes %>%
  left_join(y = Methyl_Cluster_info, by = join_by(ensembl == Ensembl)) %>%
  subset(Reg_element == 'Enhancer' & Cluster == 'C4' & GENEBIOTYPE == 'protein_coding')
rm(Methyl_Cluster_info)


##### With PReCIS-seq data #####
PreciseSeqFC_table <- read.xlsx("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/Precise2024_FoldChangeTable_v2.xlsx")
PreciseSeqFC_table$Category <- factor(PreciseSeqFC_table$Category,
                                      levels = c('upPI_nsGB', 'nsPI_upGB', 'downPI_upGB', 'downPI_nsGB', 'others'))

Peak_distance <- read.table("DataFiles/Methyl_Annotation/ClosestEpic2Peak_toPromoters.bed", sep = "\t") %>%
  select(c(V4, V10, V11, V13)) %>% rename(GeneID = V4, PeakScore = V11, Distance = V13, PeakID = V10) %>% filter(PeakScore > 0) %>%
  merge(PreciseSeqFC_table, ., by = 'GeneID', all = FALSE) %>%
  mutate(LogDistance = case_when(Distance > 0 ~ log10(Distance),
                                 Distance == 0 ~ 0,
                                 Distance < 0 ~ -log10(-Distance))) %>%
  mutate(Reg_element = if_else(abs(Distance) <= 500, 'Promoter', 'Enhancer'))

Reg_genes <- subset(Peak_distance, Category != 'others')

## Double Check before Overwrite!
write.xlsx(Reg_genes, file = 'DataFiles/GeneAnalysis_Promoter_Enhancer/PReCIS_Reg_genes.xlsx', overwrite = FALSE)

#### Check promoter vs enhancer regulated
Reg_genes <- read.xlsx('DataFiles/GeneAnalysis_Promoter_Enhancer/PReCIS_Reg_genes.xlsx')

promoter_reg_genes <- subset(Reg_genes, Reg_element == 'Promoter' & GENEBIOTYPE == 'protein_coding')
promoter_reg_genes$Symbol

dReg_regions <- read.table('/Users/kb744/Tumbar Lab Local/Precise-seq_2024/dREG/KO_exclusive_dREGs.bed')
names(dReg_regions) <- c('seqnames', 'start', 'end', 'name', 'score', 'strand', 'pvalue', 'center')
dReg_regions <- GRanges(dReg_regions)

all_genes$Reg_bydREG <- countOverlaps(all_promoters_10k, dReg_regions) > 0
dREG_regulation <- as.data.frame(all_genes) %>% select(gene_id, Reg_bydREG)

Reg_genes <- left_join(Reg_genes, dREG_regulation, by = join_by(GeneID == gene_id))
enhancer_reg_genes <- subset(Reg_genes, Reg_bydREG & GENEBIOTYPE == 'protein_coding' & Reg_element == 'Enhancer')
enhancer_reg_genes$Symbol

go <- clusterProfiler::enrichGO(enhancer_reg_genes$GeneID, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "ENSEMBL", pvalueCutoff = 0.05, readable = TRUE)
go <- as.data.frame(go)

write.xlsx(go, "DataFiles/GOAnalyses/GO_dREG_regulated_genes.xlsx")

rm(dReg_regions, dREG_regulation, Reg_genes, promoter_reg_genes, enhancer_reg_genes, go)

#### Plot GO as barplot
go_table <- read.xlsx("DataFiles/GOAnalyses/GO_dREG_regulated_genes.xlsx")

selectedGO <- c('GO:0019058', 'GO:0006457', 'GO:0019079', 'GO:0050657', 'GO:0050658')

df_to_plot <- go_table %>% arrange(Count) %>% mutate(Description = factor(Description, levels = Description)) %>%
  mutate(LogPadjust = -log10(p.adjust)) %>%
  subset(ID %in% selectedGO)

ggplot(df_to_plot, aes(x = Description, y = Count, fill = LogPadjust)) +
  geom_bar(stat = 'identity', color = 'black') +
  coord_flip() +
  scale_fill_gradient(low = 'white', high = 'darkgreen', limits = c(0, 2)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(x = NULL, y = "Gene count", fill = "-log10(p.adjust)", title = 'dREG-regulated genes') +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 12, color = 'black', face = 'plain'))


