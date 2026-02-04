library(tidyverse)
library(tidyr)
library(Seurat)
library(patchwork)
library(openxlsx)
library(clusterProfiler)
library(org.Mm.eg.db)
library(simplifyEnrichment)


## DE analysis for each cluster
CB_GO_Analyze <- function(gList, name, wb){
  require(clusterProfiler)
  require(org.Mm.eg.db)
  require(openxlsx)
  
  cat(paste0("Analyzing cluster: ", name, "\n"))
  
  genes_up <- unique(subset(gList, p_val_adj < 0.05 & avg_log2FC > 1, select = GeneID, drop = TRUE))
  genes_down <- unique(subset(gList, p_val_adj < 0.05 & avg_log2FC < -1, select = GeneID, drop = TRUE))
  
  go_up <- enrichGO(genes_up, universe = all.features, OrgDb = org.Mm.eg.db,
                    ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
  go_down <- enrichGO(genes_down, universe = all.features, OrgDb = org.Mm.eg.db,
                      ont = "ALL", keyType = "SYMBOL", pvalueCutoff = 0.05, readable = TRUE)
  
  addWorksheet(wb, paste0(name, "_Up"))
  addWorksheet(wb, paste0(name, "_Down"))
  writeData(wb, sheet = paste0(name, "_Up"), as.data.frame(go_up), rowNames = TRUE)
  writeData(wb, sheet = paste0(name, "_Down"), as.data.frame(go_down), rowNames = TRUE)
}

wb <- createWorkbook()

CB_GO_Analyze(ClusterDEgenes[['0']], name = 'Basal_1', wb)
CB_GO_Analyze(ClusterDEgenes[['2']], name = 'Basal_2', wb)
CB_GO_Analyze(ClusterDEgenes[['5']], name = 'PrePlacode', wb)
CB_GO_Analyze(ClusterDEgenes[['3']], name = 'HairFollicle', wb)
CB_GO_Analyze(ClusterDEgenes[['1']], name = 'Diff_1', wb)
CB_GO_Analyze(ClusterDEgenes[['6']], name = 'Diff_2', wb)
CB_GO_Analyze(ClusterDEgenes[['7']], name = 'Diff_3', wb)

saveWorkbook(wb, file = "ExcelFiles/GO_DEgenes_v3.xlsx", overwrite = TRUE)
rm(wb)

## Select GO terms for plotting
CB_plot_GO <- function(GO_terms, go_table){
  df_to_plot <- go_table %>% arrange(Count) %>% mutate(Description = factor(Description, levels = Description)) %>%
    mutate(LogPadjust = -log10(p.adjust)) %>%
    subset(ID %in% GO_terms)
  ggplot(df_to_plot, aes(x = Description, y = Count, fill = LogPadjust)) +
    geom_bar(stat = 'identity', color = 'black') +
    coord_flip() +
    scale_fill_gradient(low = 'white', high = 'darkgreen', limits = c(0, 8)) +
    scale_y_continuous(expand = c(0.01, 0)) +
    labs(x = NULL, y = "Gene count", fill = "-log10(p.adjust)") +
    theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
          axis.line = element_blank(),
          axis.text.y = element_text(size = 12, color = 'black', face = 'plain'),
          title = element_text(size = 10),
          aspect.ratio = length(GO_terms)/4)
}

CB_Temp_Function <- function(GO_terms, GO_table){
  
  p1 <- CB_plot_GO(GO_terms = GO_terms[['Up']], go_table = GO_table[['Up']]) +
    labs(title = "Upregulated") +
    theme(axis.title.x = element_blank(),
          legend.position = 'none')
  
  p2 <- CB_plot_GO(GO_terms = GO_terms[['Down']], go_table = GO_table[['Down']]) +
    labs(title = "Downregulated")
  
  return(p1/p2)
}

# Basal_1
GO_terms_to_plot <- list('Up' = c("GO:0019318", "GO:0140014"),
                         'Down' = c("GO:0030198", "GO:0048246", "GO:0030216", "GO:0010810"))
go_table <- list('Up' = read.xlsx("ExcelFiles/GO_DEgenes_v3.xlsx", sheet = 'Basal_1_Up', rowNames = TRUE),
                 'Down' = read.xlsx("ExcelFiles/GO_DEgenes_v3.xlsx", sheet = 'Basal_1_Down', rowNames = TRUE))
CB_Temp_Function(GO_terms_to_plot, go_table)

# Basal_2
GO_terms_to_plot <- list('Up' = c("GO:0019318", "GO:0001666"),
                         'Down' = c("GO:0030198"))
go_table <- list('Up' = read.xlsx("ExcelFiles/GO_DEgenes_v3.xlsx", sheet = 'Basal_2_Up', rowNames = TRUE),
                 'Down' = read.xlsx("ExcelFiles/GO_DEgenes_v3.xlsx", sheet = 'Basal_2_Down', rowNames = TRUE))
CB_Temp_Function(GO_terms_to_plot, go_table)



#### GO cluster analysis using simplifyEnrichment package
library(simplifyEnrichment)

go_table <- list('Up' = read.xlsx("ExcelFiles/GO_DEgenes_v3.xlsx", sheet = 'Basal_2_Up', rowNames = TRUE),
                 'Down' = read.xlsx("ExcelFiles/GO_DEgenes_v3.xlsx", sheet = 'Basal_2_Down', rowNames = TRUE))
local({
  ids <- subset(go_table[['Up']], ONTOLOGY == 'BP', select = ID, drop = TRUE)
  mat = GO_similarity(ids)
  simplifyGO(mat)
})

go_table <- list()
reg <- 'Up'
for (name in c('Basal_1', 'Basal_2', 'Basal_3', 'HairFollicle', 'Diff_1', 'Diff_2', 'Diff_3')){
  go_table[[name]] <- read.xlsx("ExcelFiles/GO_DEgenes_v3.xlsx", sheet = paste0(name, "_", reg), rowNames = TRUE)
}

go_table <- lapply(go_table, function(dt){subset(dt, ONTOLOGY =='BP')})

pdf("Plots/GO_Clusters_UpGenes.pdf", width = 9, height = 6)
simplifyGOFromMultipleLists(lt = go_table[-4], padj_cutoff = 0.05)
dev.off()

go_table <- list()
reg <- 'Down'
for (name in c('Basal_1', 'Basal_2', 'Basal_3', 'HairFollicle', 'Diff_1', 'Diff_2', 'Diff_3')){
  go_table[[name]] <- read.xlsx("ExcelFiles/GO_DEgenes_v3.xlsx", sheet = paste0(name, "_", reg), rowNames = TRUE)
}

go_table <- lapply(go_table, function(dt){subset(dt, ONTOLOGY =='BP')})

pdf("Plots/GO_Clusters_DownGenes.pdf", width = 9, height = 6)
simplifyGOFromMultipleLists(lt = go_table[-4], padj_cutoff = 0.05)
dev.off()


#### Make Dotplot for scRNA-seq DE genes in normal populations ####
wb = loadWorkbook("ExcelFiles/GO_DEgenes_v3.xlsx")
GO_tables <- list()
for (name in sheets(wb)){
  GO_tables[[name]] <- read.xlsx(wb, name, rowNames = TRUE)
}
GO_tables[['HairFollicle_Up']] <- NULL
GO_tables[['HairFollicle_Down']] <- NULL
rm(wb, name)

### Manually select GO terms from GO_DEgenes_v3.xlsx for plotting
GO_terms_to_plot <- c("GO:0006007", "GO:0006735", "GO:0098813", "GO:0006096", "GO:0048285", "GO:0007051", "GO:0044839", "GO:0051321", #Basal_1_Up
                      "GO:0071456",
                      "GO:0009060", "GO:0006119", "GO:0022900", "GO:0140241", "GO:0008380", "GO:0006417",
                      "GO:0005200", # Diff_1_Up
                      "GO:0098984",
                      "GO:0030198", "GO:0030216", "GO:0010810", # Basal_1_Down
                      "GO:0062023",
                      "GO:0008356",
                      "GO:0005604", # Diff_1_Down
                      "GO:0035113", "GO:0032330", "GO:0001708", "GO:0022405", "GO:0048730",
                      "GO:0003334"
                      )

go_table <- bind_rows(GO_tables, .id = 'GeneList') %>% subset(ID %in% GO_terms_to_plot) %>%
  mutate(ID = factor(ID, levels = GO_terms_to_plot)) %>% arrange(desc(ID)) %>%
  mutate(Description = factor(Description, levels = unique(Description))) %>%
  mutate(GeneList = factor(GeneList, levels = c("Basal_1_Up", "Basal_2_Up", "Basal_3_Up", 'Diff_1_Up', "Diff_2_Up", "Diff_3_Up",
                                                "Basal_1_Down", "Basal_2_Down", "Basal_3_Down", 'Diff_1_Down', "Diff_2_Down", "Diff_3_Down")))

ggplot(go_table, aes(x = GeneList, y = Description, color = qvalue)) +
  geom_point(aes(size = Count)) +
  scale_size(range = c(3,9), breaks = c(4, 7, 10)) +
  scale_color_gradient(low = 'darkgreen', high = 'gray90', limits = c(0, 0.05)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_text(color = 'black', face = 'plain'),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

