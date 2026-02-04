library(Seurat)
library(tidyverse)


harmonyEpiCells <- readRDS(file = "Harmony_integrated_epidermal_w_Merkel.RData")

gene_count_summary <- table(
  unlist(
    lapply(ClusterDEgenes, function(df){
      lt <- df$GeneID
    })
  )
)
gene_count_summary <- as.data.frame(gene_count_summary)

ggplot(gene_count_summary, aes(x = Freq)) +
  geom_histogram()
sum(gene_count_summary$Freq > 6)

Idents(harmonyEpiCells) <- 'Identity'

cluster_cnt <- table(Idents(harmonyEpiCells))
gene_cnt <- lapply(ClusterDEgenes, function(df){
  df <- subset(df, avg_log2FC > 1 & p_val_adj < 0.05)
  cnt <- data.frame(GeneID = df$GeneID)
})

names(gene_cnt) <- Cluster_idents[names(gene_cnt)]

for (name in names(gene_cnt)) {
  gene_cnt[[name]]$Freq <- cluster_cnt[[name]]
  colnames(gene_cnt[[name]])[2] <- name
}

gene_count_summary <- reduce(gene_cnt, full_join, by = 'GeneID') %>% subset(!duplicated(GeneID))
gene_count_summary[is.na(gene_count_summary)] <- 0
gene_count_summary$Sum <- rowSums(gene_count_summary[,2:11])
ggplot(gene_count_summary, aes(x = Sum)) +
  geom_histogram() +
  labs(x = '#Cells genes are upregulated')
sum(gene_count_summary$Sum > 6000)



