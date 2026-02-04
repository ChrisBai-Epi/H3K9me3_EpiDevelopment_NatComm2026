library(openxlsx)
library(tidyverse)

### This scirpt uses data table generated from ProfilingMethyAroundGenes_2.R
methylInfo_promoters <- read.xlsx('DataFiles/Methyl_Annotation/Methyl_annotation_allRegions_allGenes.xlsx', sheet = 'promoter')
methylInfo_flank10k <- read.xlsx('DataFiles/Methyl_Annotation/Methyl_annotation_allRegions_allGenes.xlsx', sheet = 'flank10k')

## For nascent up genes vs random control
Nascent_up_genes
random_Genes

methylInfo <- merge(methylInfo_promoters[, c('ensembl', 'name', 'MeanIntensity')],
                    methylInfo_flank10k[, c('ensembl', 'name', 'MeanIntensity')],
                    by = c('ensembl', 'name'))
df_to_plot <- methylInfo %>% rename(promoter = MeanIntensity.x, flank = MeanIntensity.y) %>%
  mutate(Group = case_when(ensembl %in% Nascent_up_genes ~ 'Nascent_up',
                           ensembl %in% random_Genes ~ 'Random',
                           .default = 'others')) %>%
  subset(Group != 'others') %>%
  pivot_longer(cols = c('promoter', 'flank'), names_to = 'Region', values_to = 'MeanIntensity') %>%
  mutate(Region = factor(Region, levels = c('promoter', 'flank')))

ggplot(df_to_plot, aes(x = Region, y = MeanIntensity, fill = Group)) +
  geom_boxplot(outliers = FALSE, linewidth = 1, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c('darkgreen', '#B6C69B')) +
  labs(x = NULL)

ggplot(df_to_plot, aes(x = Region, y = MeanIntensity, fill = Group)) +
  geom_bar(stat = 'summary', fun = 'mean', color = 'black', size = 1, position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, size = 1, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c('darkgreen', '#B6C69B')) +
  labs(x = NULL)

values_nascent_up <- subset(methylInfo, ensembl %in% Nascent_up_genes)
values_random <- subset(methylInfo, ensembl %in% random_Genes)

t.test(values_nascent_up$MeanIntensity.x, values_random$MeanIntensity.x)
t.test(values_nascent_up$MeanIntensity.y, values_random$MeanIntensity.y)

rm(values_nascent_up, values_random)

## For the four category of genes
names(Nascent_GOI)

local({
  set.seed(6308)
  rdm_df <- sample_n(as.data.frame(all_genes), 150)
  random_Genes_150 <<- rdm_df$gene_id
})

df_to_plot <- methylInfo %>% rename(promoter = MeanIntensity.x, flank = MeanIntensity.y) %>%
  mutate(Group = case_when(ensembl %in% Nascent_GOI$upPI_nsGB ~ 'upPI_nsGB',
                           ensembl %in% Nascent_GOI$nsPI_upGB ~ 'nsPI_upGB',
                           ensembl %in% Nascent_GOI$downPI_upGB ~ 'downPI_upGB',
                           ensembl %in% Nascent_GOI$downPI_nsGB ~ 'downPI_nsGB',
                           ensembl %in% methylated_TSS$gene_id ~ 'methylated_TSS',
                           ensembl %in% random_Genes_150 ~ 'random',
                           .default = 'others')) %>%
  mutate(Group = factor(Group, c('upPI_nsGB', 'nsPI_upGB', 'downPI_upGB', 'downPI_nsGB', 'random', 'methylated_TSS', 'others'))) %>%
  subset(Group != 'others')

ggplot(df_to_plot, aes(x = Group, y = promoter)) +
  geom_boxplot(outliers = FALSE, linewidth = 1, fill = 'darkgreen') +
  labs(x = NULL, y = 'Methylation signal')

ggplot(df_to_plot, aes(x = Group, y = promoter)) +
  geom_bar(stat = 'summary', fun = 'mean', color = 'black', fill = 'darkgreen', size = 1, position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, size = 1, position = position_dodge(width = 0.8)) +
  coord_cartesian(ylim = c(-0.055, 0.08)) +
  labs(x = NULL, y = 'Methylation signal')

values_upPInsGB <- subset(methylInfo, ensembl %in% Nascent_GOI$upPI_nsGB)
values_downPInsGB <- subset(methylInfo, ensembl %in% Nascent_GOI$downPI_nsGB)
values_random <- subset(methylInfo, ensembl %in% random_Genes_150)
values_methylated <- subset(methylInfo, ensembl %in% methylated_TSS$gene_id)

t.test(values_upPInsGB$MeanIntensity.x, values_random$MeanIntensity.x)
t.test(values_downPInsGB$MeanIntensity.x, values_random$MeanIntensity.x)
t.test(values_upPInsGB$MeanIntensity.x, values_methylated$MeanIntensity.x)

rm(values_upPInsGB, values_random, values_downPInsGB, values_methylated)

## Barplot for comparing 3 nascent categories
names(Nascent_GOI)

local({
  set.seed(6308)
  rdm_df <- sample_n(as.data.frame(all_genes), 150)
  random_Genes_150 <<- rdm_df$gene_id
})

df_to_plot <- methylInfo %>% rename(promoter = MeanIntensity.x, flank = MeanIntensity.y) %>%
  mutate(Group = case_when(ensembl %in% Nascent_GOI$upPI_nsGB ~ 'Category_1',
                           ensembl %in% Nascent_GOI$nsPI_upGB ~ 'Category_2',
                           ensembl %in% Nascent_GOI$downPI_upGB | ensembl %in% Nascent_GOI$downPI_nsGB ~ 'Category_3',
                           ensembl %in% random_Genes_150 ~ 'Random',
                           .default = 'Others')) %>%
  mutate(Group = factor(Group, c('Category_1', 'Category_2', 'Category_3', 'Random', 'Others'))) %>%
  subset(Group != 'Others')

ggplot(df_to_plot, aes(x = Group, y = promoter)) +
  geom_bar(stat = 'summary', fun = 'mean', color = 'black', fill = 'darkgreen', size = 1, position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, width = 0.2, size = 1, position = position_dodge(width = 0.8)) +
  coord_cartesian(ylim = c(-0.055, 0)) +
  labs(x = NULL, y = 'Methylation signal')

values_Category_1 <- subset(methylInfo, ensembl %in% Nascent_GOI$upPI_nsGB)
values_Category_3 <- subset(methylInfo, ensembl %in% Nascent_GOI$downPI_nsGB | ensembl %in% Nascent_GOI$downPI_upGB)
values_random <- subset(methylInfo, ensembl %in% random_Genes_150)

t.test(values_Category_1$MeanIntensity.x, values_random$MeanIntensity.x) ## pvalue=0.5376
t.test(values_Category_3$MeanIntensity.x, values_random$MeanIntensity.x) ## pvalue=0.0912

rm(values_Category_1, values_random, values_Category_3)

