library(BRGenomics)
library(GenomicRanges)
library(tidyverse)

#Import Bam files
samples_list <- import_bam_PROseq(
  list(
    "../BamFiles/CT_Input_A_deduplicated.bam",
    "../BamFiles/CT_Input_B_deduplicated.bam",
    "../BamFiles/CT_IP_A_deduplicated.bam",
    "../BamFiles/CT_IP_B_deduplicated.bam",
    "../BamFiles/KO_Input_A_deduplicated.bam",
    "../BamFiles/KO_Input_B_deduplicated.bam",
    "../BamFiles/KO_IP_A_deduplicated.bam",
    "../BamFiles/KO_IP_B_deduplicated.bam"
  ), mapq = 0, ncores = 20,
  trim.to = '5p'
)

names(samples_list) <- c("CT_Input_A", "CT_Input_B", "CT_IP_A", "CT_IP_B",
                         "KO_Input_A", "KO_Input_B", "KO_IP_A", "KO_IP_B")

#Read Annotation File
Annotations <- read.table("RefSeq_MostActiveTranscriptPerGene_110123.bed", 
                          header=FALSE, stringsAsFactors = FALSE)
colnames(Annotations) <- c("chrom", "start", "end", "width", "strand", "transcript_id", "gene_name", "gene_type", "gene_id")

#Optionally filter for genes (for example, protein-coding genes that are longer than 2.5kb)
CodingGenes_filtered <- subset(Annotations, subset = gene_type == "protein_coding" & width > 2500)
Genes_filtered <- subset(Annotations, subset = width > 2500)

Genes_Granges <- sort(GRanges(Genes_filtered))


#Read count table by STAR
StarCounts_list <- list()
for (file in list.files("../StarCounts", full.names = TRUE)){
  name <- sub(".*/(.*)_ReadsPerGene\\.out\\.tab", "\\1", file)
  StarCounts_list[[name]] <- read.table(file = file, skip = 4) %>% dplyr::select(V1, V4)
  names(StarCounts_list[[name]]) <- c("ID", name)
}
rm(file, name)

StarCounts <- Reduce(function(x, y) merge(x , y, by = 'ID', all = TRUE), StarCounts_list)
StarCounts$ID <- sub("\\..*", "", StarCounts$ID)

StarCounts_Fly <- StarCounts[!grepl("^ENSMUSG", StarCounts$ID),]
StarCounts_Mouse <- StarCounts[grep("^ENSMUSG", StarCounts$ID),]


#Define normalizing factors
readCountSummary <- getSpikeInCounts(samples_list, si_pattern = '^dm_', ncores = 10)

NF_RPM <- 1e6 / readCountSummary$exp_reads
NF_SpikeIn <- 1e3 / readCountSummary$spike_reads

###Calculate scaling factors for generating spike-in normalized bigwig files
write.csv(readCountSummary, file = "ReadCountSummary.csv")

#Get counts or promoter proximal region and gene body regions
PauseCounts <- getCountsByRegions(samples_list, genebodies(Genes_Granges, 0, 250, fix.end = "start")) %>%
  add_column(GeneID = Genes_Granges$gene_id, .before = 1) %>% column_to_rownames("GeneID")

GeneBodyCounts <- getCountsByRegions(samples_list, genebodies(Genes_Granges, 500, -500)) %>%
  add_column(GeneID = Genes_Granges$gene_id, .before = 1) %>% column_to_rownames("GeneID")

FullCounts <- getCountsByRegions(samples_list, Genes_Granges) %>%
  add_column(GeneID = Genes_Granges$gene_id, .before = 1) %>% column_to_rownames("GeneID")

#Make sure PauseCounts and GeneBodyCounts have the same number of rows and row order
all(rownames(GeneBodyCounts) == rownames(PauseCounts))

#######Normalization######################
gene_length <- width(genebodies(Genes_Granges, 500, -500)) / 1e3
GeneBodyCounts_RPKM <- GeneBodyCounts * outer(1/gene_length, NF_RPM)
GeneBodyCounts_SPKK <- GeneBodyCounts * outer(1/gene_length, NF_SpikeIn)

PauseCounts_RPKM <- PauseCounts * outer(1000/250, NF_RPM)
PauseCounts_SPKK <- PauseCounts * outer(1000/250, NF_SpikeIn)

gene_full_length <- width(Genes_Granges) / 1e3
FullCounts_SPKK <- FullCounts * outer(1/gene_full_length, NF_SpikeIn)


###Plot to check if there is a global change
library(ggpubr)

df_to_plot <- GeneBodyCounts_RPKM %>%  select(c(3,4,7,8)) %>% filter(rowSums(Merged_Count_Table) >= 80) %>% rownames_to_column(var = 'ID') %>%
  pivot_longer(cols = 2:5, names_to = 'Sample', values_to = 'RPKM') %>%
  mutate(Genotype = case_when(grepl("^CT", Sample) ~ 'Ctrl', grepl("^KO", Sample) ~ 'TKO'))
y.title <- 'RPKM'

df_to_plot <- GeneBodyCounts_SPKK %>%  select(c(3,4,7,8)) %>% filter(rowSums(Merged_Count_Table) >= 80) %>% rownames_to_column(var = 'ID') %>%
  pivot_longer(cols = 2:5, names_to = 'Sample', values_to = 'RPKM') %>%
  mutate(Genotype = case_when(grepl("^CT", Sample) ~ 'Ctrl', grepl("^KO", Sample) ~ 'TKO'))
y.title <- 'Spike-in Norm. Counts'

ggplot(df_to_plot, aes(x = Sample, y = RPKM, fill = Genotype)) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = genoColor) +
  labs(x = element_blank(),
       y = y.title,
       fill = 'Sample')

##Using merged table
df_to_plot <- normalizedCounts_RPKM %>% filter(rowSums(Merged_Count_Table) >= 80) %>% rownames_to_column(var = 'ID') %>%
  pivot_longer(cols = 2:3, names_to = 'Sample', values_to = 'RPKM') %>%
  mutate(Genotype = case_when(grepl("^CT", Sample) ~ 'Ctrl', grepl("^KO", Sample) ~ 'TKO'))
y.title <- 'RPKM'

df_to_plot <- normalizedCounts_SPKK %>% filter(rowSums(Merged_Count_Table) >= 80) %>% rownames_to_column(var = 'ID') %>%
  pivot_longer(cols = 2:3, names_to = 'Sample', values_to = 'RPKM') %>%
  mutate(Genotype = case_when(grepl("^CT", Sample) ~ 'Ctrl', grepl("^KO", Sample) ~ 'TKO'))
y.title <- 'Spike-in Norm. Counts'

ggplot(df_to_plot, aes(x = Genotype, y = RPKM, fill = Genotype)) +
  geom_boxplot(outliers = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = genoColor) +
  labs(x = element_blank(),
       y = y.title)

t.spikeIn <- t.test(df_to_plot[df_to_plot$Genotype == 'Ctrl', 3], df_to_plot[df_to_plot$Genotype == 'TKO', 3])
t.RPKM <- t.test(df_to_plot[df_to_plot$Genotype == 'Ctrl', 3], df_to_plot[df_to_plot$Genotype == 'TKO', 3])

#T test for the spike-in normalized total reads
Reads_per_spikeIn <- readCountSummary$exp_reads / readCountSummary$spike_reads

df_summary <- data.frame(Sample = rep(c("Ctrl", "TKO"), each = 2), RPS = Reads_per_spikeIn[c(3,4,7,8)]) %>%
  group_by(Sample) %>% summarise(Mean = mean(RPS), SE = sd(RPS) / sqrt(n()))

ggplot(df_summary, aes(x = Sample, y = Mean, fill = Sample)) +
  geom_bar(stat = 'identity', width = 0.8, show.legend = FALSE) + 
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  scale_fill_manual(values = genoColor) +
  labs(x = element_blank(), y = "Read counts per spike-in read") +
  theme(axis.text.x = element_text(size = 14))

t.total <- t.test(Reads_per_spikeIn[3:4], Reads_per_spikeIn[7:8])

#Using merged normalization values
t.spikeIn <- t.test(normalizedCounts_SPKK$CT_IP_SPKK, normalizedCounts_SPKK$KO_IP_SPKK)
t.RPKM <- t.test(normalizedCounts_RPKM$CT_IP_RPKM, normalizedCounts_RPKM$KO_IP_RPKM)

###Generate base-pair resolution bigwig files
mm_chrs <- c(paste0("chr", 1:19), "chrX", "chrY", "chrM")

Merge_and_export_bw <- function(gList, baseName) {
  require(rtracklayer)
  gr <- mergeGRangesData(gList)
  gr <- subset(gr, seqnames(gr) %in% mm_chrs)
  seqlevels(gr) <- mm_chrs
  gr_fwd <- subset(gr, strand == '+')
  gr_rev <- subset(gr, strand == '-')
  score(gr_rev) <- -score(gr_rev)
  
  export.bw(gr_fwd, con = paste0(baseName, "_fwd.bw"))
  export.bw(gr_rev, con = paste0(baseName, "_rev.bw"))
}

Merge_and_export_bw(samples_list[1:2], baseName = '../BrowserTracks/BasePairResoBigWigs/CT_Input')
Merge_and_export_bw(samples_list[3:4], baseName = '../BrowserTracks/BasePairResoBigWigs/CT_IP')
Merge_and_export_bw(samples_list[5:6], baseName = '../BrowserTracks/BasePairResoBigWigs/KO_Input')
Merge_and_export_bw(samples_list[7:8], baseName = '../BrowserTracks/BasePairResoBigWigs/KO_IP')


##Annotation and ggplot
file.edit("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/UtilityFunctions.R")

rm(list = ls()[grep("^gr", ls())])

save.image()
