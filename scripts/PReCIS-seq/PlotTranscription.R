library(openxlsx)
library(tidyverse)

PreciseSeqFC_table <- read.xlsx("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/Genomics_DE_files/Precise2024_FoldChangeTable_annotated.xlsx")
PreciseSeqFC_table <- PreciseSeqFC_table %>%
  mutate(Category = if_else(Reg_PI == 'down', 'DownPI_upGB', Category))

Nascent_GOI <- list()
for (name in unique(PreciseSeqFC_table$Category)){
  Nascent_GOI[[name]] <- subset(PreciseSeqFC_table, subset = Category == name, select = GeneID, drop = TRUE)
  write.table(Nascent_GOI[[name]], file = paste0('MotifAnalysisFiles/PreciseSeq_', name, '.txt'), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

for (set in names(Nascent_GOI)){
  gr <- subset(AllGenes, gene_id %in% Nascent_GOI[[set]])
  df <- as.data.frame(gr)
  df$start <- df$start - 1
  df$score <- 0
  bed_df <- df[, c('seqnames', 'start', 'end', 'gene_id', 'score', 'strand', 'gene_name', 'gene_type')]
  write.table(bed_df, file = paste0("ProfilePlots/Region_", set, ".bed"), quote = FALSE, sep = '\t', col.names = FALSE, row.names = FALSE)
}
rm(set, gr, df, bed_df)

command <- paste0("computeMatrix scale-regions -S ",
                  "SpikeIn_Norm_BigWigs/CT_IP_fwd.bw SpikeIn_Norm_BigWigs/CT_IP_rev.bw ",
                  "SpikeIn_Norm_BigWigs/KO_IP_fwd.bw SpikeIn_Norm_BigWigs/KO_IP_rev.bw ",
                  "-R ProfilePlots/Region_Up_PI_ns_GB_genes.bed ProfilePlots/Region_Ns_PI_Up_GB_genes.bed ",
                  "ProfilePlots/Region_Down_PI_Up_GB_genes.bed ProfilePlots/Region_Down_PI_ns_GB_genes.bed ",
                  "ProfilePlots/Region_Others.bed ",
                  "-o ProfilePlots/Matrix_Nascent.gz ", 
                  "--regionBodyLength 4000 --binSize 20 --sortRegions descend ",
                  "-p 10 -b 2000 -a 2000")
system(command)
command <- paste0("plotProfile -m ProfilePlots/Matrix_Nascent.gz ", 
                  "--samplesLabel CT-fwd CT-rev KO-fwd KO-rev ",
                  "-o ProfilePlots/ProfilePlot_nascentRNA.pdf ",
                  "--regionsLabel upPI_nsGB nsPI_upGB downPI_upGB downPI_nsGB others ",
                  "--yMin -5 -5 -10 -5 -5 --yMax 5 5 10 5 5 ",
                  "--perGroup --averageType mean --numPlotsPerRow 1 --plotHeight 6 --plotWidth 9 --legendLocation upper-center")
system(command)


### Plot bargraph for GB counts, PP counts and Pausing Index

df_to_plot <- FC_table_w_test %>% select(c(Symbol, Reg_PI, Reg_GB)) %>% 
  cbind(., Merged_PI_Table[row.names(.), 1:6]) %>% rownames_to_column(var = "Ensembl") %>% select(!Symbol)
colnames(df_to_plot)[4:9] <- c("Ctrl_GB", "TKO_GB", "Ctrl_PP", "TKO_PP", "Ctrl_PI", "TKO_PI")
df_to_plot <- df_to_plot %>%
  mutate(Category = paste0(Reg_PI, "PI_", Reg_GB, "GB")) %>% 
  mutate(Category = if_else(Category %in% c("nsPI_nsGB", "upPI_upGB", "nsPI_downGB"), "others", Category)) %>%
  pivot_longer(cols = 4:9, names_to = c("Genotype", "DataType"), names_sep = "_", values_to = "Value") %>%
  pivot_wider(names_from = "DataType", values_from = "Value")

df_to_plot$Category <- factor(df_to_plot$Category, levels = c('upPI_nsGB', 'nsPI_upGB', 'downPI_upGB', 'downPI_nsGB', 'others'))


p1 <- ggplot(df_to_plot, aes(x = Genotype, y = GB, fill = Genotype)) +
  geom_bar(stat = 'summary', fun = 'mean', color = 'black', size = 1, show.legend = FALSE) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, linewidth = 1) +
  scale_fill_manual(values = genoColor) +
  facet_wrap(~Category, ncol = 1, scales = 'free_y', strip.position = 'left') +
  labs(x = NULL) +
  theme(strip.text = element_text(), panel.spacing = unit(1, 'cm'), axis.text.x = element_text())

p2 <- ggplot(df_to_plot, aes(x = Genotype, y = PP, fill = Genotype)) +
  geom_bar(stat = 'summary', fun = 'mean', color = 'black', size = 1, show.legend = FALSE) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, linewidth = 1) +
  scale_fill_manual(values = genoColor) +
  facet_wrap(~Category, ncol = 5) +
  labs(x = NULL) +
  theme(strip.text = element_blank(), panel.spacing = unit(1, 'cm'), axis.text.x = element_blank())

p3 <- ggplot(df_to_plot, aes(x = Genotype, y = PI, fill = Genotype)) +
  geom_bar(stat = 'summary', fun = 'mean', color = 'black', size = 1, show.legend = FALSE) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, linewidth = 1) +
  scale_fill_manual(values = genoColor) +
  facet_wrap(~Category, ncol = 5) +
  labs(x = NULL) +
  theme(strip.text = element_blank(), panel.spacing = unit(1, 'cm'), 
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5))

p1 / p2 / p3


#### Plot: categories as columns, regions as rows, sharing the same y scale among categroies
df_to_plot <- FC_table_w_test %>% select(c(Symbol, Reg_PI, Reg_GB)) %>% 
  cbind(., Merged_PI_Table[row.names(.), 1:6]) %>% rownames_to_column(var = "Ensembl") %>% select(!Symbol)
colnames(df_to_plot)[4:9] <- c("Ctrl_GB", "TKO_GB", "Ctrl_PP", "TKO_PP", "Ctrl_PI", "TKO_PI")
df_to_plot <- df_to_plot %>%
  mutate(Category = paste0(Reg_PI, "PI_", Reg_GB, "GB")) %>% 
  mutate(Category = if_else(Category %in% c("nsPI_nsGB", "upPI_upGB", "nsPI_downGB"), "others", Category)) %>%
  pivot_longer(cols = 4:9, names_to = c("Genotype", "DataType"), names_sep = "_", values_to = "Value")

df_to_plot$Category <- factor(df_to_plot$Category, levels = c('upPI_nsGB', 'nsPI_upGB', 'downPI_upGB', 'downPI_nsGB', 'others'))
df_to_plot$DataType <- factor(df_to_plot$DataType, levels = c('PP', 'GB', 'PI'))

ggplot(df_to_plot, aes(x = Genotype, y = Value, fill = Genotype)) +
  geom_bar(stat = 'summary', fun = 'mean', color = 'black', size = 1, show.legend = FALSE) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, linewidth = 1) +
  scale_fill_manual(values = genoColor) +
  facet_wrap(Category ~ DataType, scales = 'free_y', ncol = 3, axis.labels = 'margins') +
  labs(x = NULL, y = NULL) +
  theme(strip.text = element_blank(), panel.spacing = unit(0.5, 'cm'),
        axis.text.x = element_text())

ggplot(df_to_plot, aes(x = Genotype, y = Value, fill = Genotype)) +
  geom_bar(stat = 'summary', fun = 'mean', color = 'black', size = 1, show.legend = FALSE) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, linewidth = 1) +
  scale_fill_manual(values = genoColor) +
  facet_grid(rows = vars(DataType), cols = vars(Category), scales = 'free_y') +
  labs(x = NULL, y = NULL) +
  theme(panel.spacing = unit(0.5, 'cm'),
        axis.text.x = element_text())

#### Plot Box and whisker plot
ggplot(df_to_plot, aes(x = Genotype, y = Value, fill = Genotype)) +
  geom_boxplot(outliers = FALSE, linewidth = 0.5) +
  scale_fill_manual(values = genoColor) +
  facet_grid(rows = vars(DataType), cols = vars(Category), scales = 'free_y') +
  labs(x = NULL, y = NULL) +
  theme(panel.spacing = unit(0.5, 'cm'),
        axis.text.x = element_text())


#### Plot: same barplot as above for when downPI_upGB/downPI_nsGB are merged
df_to_plot <- FC_table_w_test %>% select(c(Symbol, Reg_PI, Reg_GB, Category)) %>% 
  cbind(., Merged_PI_Table[row.names(.), 1:6]) %>% rownames_to_column(var = "Ensembl") %>% select(!Symbol)
colnames(df_to_plot)[5:10] <- c("Ctrl_GB", "TKO_GB", "Ctrl_PP", "TKO_PP", "Ctrl_PI", "TKO_PI")
df_to_plot <- df_to_plot %>%
  mutate(Category = case_when(Category == 'upPI_nsGB' ~ 'Category_1',
                              Category == 'nsPI_upGB' ~ 'Category_2',
                              Category %in% c('downPI_upGB', 'downPI_nsGB') ~ 'Category_3',
                              .default = 'Others')) %>%
  pivot_longer(cols = 5:10, names_to = c("Genotype", "DataType"), names_sep = "_", values_to = "Value")

df_to_plot$Category <- factor(df_to_plot$Category, levels = c('Category_1', 'Category_2', 'Category_3', 'Others'))
df_to_plot$DataType <- factor(df_to_plot$DataType, levels = c('PP', 'GB', 'PI'))


ggplot(df_to_plot, aes(x = Genotype, y = Value, fill = Genotype)) +
  geom_bar(stat = 'summary', fun = 'mean', color = 'black', size = 1, show.legend = FALSE) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, linewidth = 1) +
  scale_fill_manual(values = genoColor) +
  facet_grid(rows = vars(DataType), cols = vars(Category), scales = 'free_y') +
  labs(x = NULL, y = NULL) +
  theme(panel.spacing = unit(0.5, 'cm'),
        axis.text.x = element_text())



