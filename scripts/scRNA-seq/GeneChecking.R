library(Seurat)
library(openxlsx)
library(org.Mm.eg.db)
library(GO.db)
library(RColorBrewer)

library(tidyverse)
library(patchwork)

harmonyEpiCells <- readRDS(file = "Harmony_integrated_epidermal_w_Merkel.RData")

### Cell cycle genes
FeaturePlot(harmonyEpiCells, features = c('Mki67', 'Top2a'), split.by = 'genotype', reduction = 'umap')
VlnPlot(harmonyEpiCells, features = c('Mki67', 'Top2a', 'Cdca3', 'Cks2', 'Ccnb2', 'Bub3'), split.by = 'genotype',
        alpha = 0)

GO_cell_cycle_genes <- AnnotationDbi::select(org.Mm.eg.db, keytype = 'GOALL', keys = 'GO:0007049', columns = 'SYMBOL')
GO_cell_cycle_genes <- unique(GO_cell_cycle_genes$SYMBOL)

harmonyEpiCells <- AddModuleScore(harmonyEpiCells, features = list(GO_cell_cycle_genes), name = "CellCycle_")

plots <- FeaturePlot(harmonyEpiCells, features = 'CellCycle_1', split.by = 'genotype')
plots[[1]]  <- plots[[1]] + scale_color_gradient2(low = 'blue', high = 'red', midpoint = 0.2)
plots[[2]]  <- plots[[2]] + scale_color_gradient2(low = 'blue', high = 'red', midpoint = 0.2)
plots[[1]] + plots[[2]]

VlnPlot(harmonyEpiCells, features = 'CellCycle_1', split.by = 'genotype')


### Antiviral genes
GO_antiviral_genes <- AnnotationDbi::select(org.Mm.eg.db, keytype = 'GOALL', keys = 'GO:0002230', columns = 'SYMBOL')
GO_antiviral_genes <- unique(GO_antiviral_genes$SYMBOL)

harmonyEpiCells <- AddModuleScore(harmonyEpiCells, features = list(GO_antiviral_genes), name = "AntiViral_")

FeaturePlot(harmonyEpiCells, features = 'AntiViral_1', split.by = 'genotype')
VlnPlot(harmonyEpiCells, features = 'AntiViral_1', split.by = 'genotype')


### Inflammation gene checking, genes from Lyu et al. Cell 2024
Antiviral_sensors <- c('Cgas', 'Sting', 'Aim2', 'Mda5', 'Rigi', 'Ddx60', 'Tlr3', 'Tlr7','Pkr', 'Ifi44')
Antiviral_effectors <- c('Irf1', 'Irf3', 'Irf7', 'Gsdma', 'Gsdma3', 'Isg15', 'Isg20', 'Oas12', 'Oas1b', 'Oas3', 'Mx1', 'Mx2')
Antiviral_induced <- c('Aim2', 'Eif2ak2', 'Ifih1', 'Ddx58', 'Dtx3l', 'Parp9', 'Zc3h12a', 'Traf3ip2', 'Creb3', 'Isg15')
Antiviral_uninduced <- c('Cgas', 'Tmem173', 'Tlr3', 'Mavs')

VlnPlot(harmonyEpiCells, features = Antiviral_sensors, split.by = 'genotype', ncol = 3)
VlnPlot(harmonyEpiCells, features = Antiviral_effectors, split.by = 'genotype', ncol = 4)
VlnPlot(harmonyEpiCells, features = Antiviral_induced, split.by = 'genotype', ncol = 3, alpha = 0)
VlnPlot(harmonyEpiCells, features = Antiviral_uninduced, split.by = 'genotype', ncol = 2, alpha = 0)

plots <- VlnPlot(harmonyEpiCells, features = Antiviral_induced, split.by = 'genotype',
                 ncol = 1, alpha = 0, combine = FALSE)
plots[1:9] <- lapply(plots[-10], function(p){
  p <- p + theme(axis.text.x = element_blank(),
                 plot.title = element_text(size = 10, hjust = 1, vjust = 0.1),
                 legend.position = 'none',
                 axis.title = element_blank(),
                 plot.margin = margin(r = 20))
})
plots[[10]] <- plots[[10]] + theme(axis.title = element_blank(),
                                 plot.title = element_text(size = 10, hjust = 1, vjust = 0.1),
                                 plot.margin = margin(r = 20))
wrap_plots(plots, ncol = 1)

## Two genes checked by immunofluorescence staining
FeaturePlot(harmonyEpiCells, features = c('Tnf', 'Ptgs2'), split.by = 'genotype', by.col = TRUE, raster = TRUE,
            pt.size = 4)

## Check expression of HMTs
VlnPlot(harmonyEpiCells, features = c('Suv39h1', 'Suv39h2', 'Setdb1'), alpha = 0, ncol = 1,
        idents = c('Basal_1', 'Basal_2', 'Basal_3', 'Diff_1', 'Diff_2', 'Diff_3', 'HairFollicle', 'Merkel'))


### Inflammation response genes from Garcia et al., Dev Cell 2024
Inflammation_factors <- c('Cxcl5', 'Cxcl10', 'Plekho', 'Cxcl14', 'Il17', 'Csf2', 'Il6', 'Inba', 'Tnfa')
DotPlot(harmonyEpiCells, features = Inflammation_factors, 
        dot.scale = 6, split.by = "genotype",
        cols = c("red", "blue"), col.min = -5, col.max = 5) +
  RotatedAxis()
VlnPlot(harmonyEpiCells, features = Inflammation_factors, split.by = 'genotype', ncol = 3, alpha = 0)


### Checking for DNA repair genes
Genes_DNA_damage <- read.xlsx("/Users/kb744/Library/CloudStorage/OneDrive-CornellUniversity/Tumbar Lab/scRNA-seq/SeuratAnalysis/ExcelFiles/DNA_Damage_ResponseGenes.xlsx")
Genes_DNA_damage$Symbol <- mapIds(org.Mm.eg.db, keys = Genes_DNA_damage$GeneID, column = 'SYMBOL', keytype = "MGI")
Genes_DNA_damage$GO_Name <- Term(GOTERM[Genes_DNA_damage$GO_ID])

GenesNames <- unique(subset(Genes_DNA_damage, GO_Name == "DNA damage response", select = "Symbol", drop = TRUE))
harmonyEpiCells <- AddModuleScore(harmonyEpiCells, features = list(GenesNames), name = "DNA_Damage_Response_")

FeaturePlot(harmonyEpiCells, features = 'DNA_Damage_Response_1', split.by = 'genotype')
VlnPlot(harmonyEpiCells, features = 'DNA_Damage_Response_1', split.by = 'genotype')

DNA_repair_genes <- c("Trp63", "Brca1", "Aatf", "Brca2", "Chaf1b", "Helb", "Chaf1a", "Exo1",
                      "Dgcr8", "Nek1", "Topbp1", "Macrod2", "Fbxo5", "Nbn", "Pole", "Vav3",
                      "Zgrf1", "Baz1b", "Spidr", "Wdr76", "Msh6", "Rad51b", "Rnf168", "Mms22l",
                      "Neil3", "Msh3", "Timeless", "Dna2", "Dtl", "Blm", "Paxip1", "Prkdc",
                      "Tank", "Rad51ap1", "Pold3", "Brip1", "Rnf138", "Abl1", "Ascc2", "Clspn",
                      "Bard1", "Slf1", "Cbx5", "Slf2", "Atad5", "Fancm", "Cradd", "Fancl",
                      "Ints3", "Eya3", "Mcrs1", "Rpa2", "Apex2", "Mre11a", "Fancc", "Rad50",
                      "Rad51", "Taok3", "Rfwd3", "Cdk2", "Mms19", "Atm", "Rad1", "Rad18")
harmonyEpiCells <- AddModuleScore(harmonyEpiCells, features = list(DNA_repair_genes), name = "DNA_Repair_")

VlnPlot(harmonyEpiCells, features = 'DNA_Repair_1', split.by = 'genotype', alpha = 0)
FeaturePlot(harmonyEpiCells, features = 'DNA_Repair_1', split.by = 'genotype', reduction = 'umap', cols = c('gray90', 'red'))


##Bona-fide Pc marker
FeaturePlot(harmonyEpiCells, features = c('Shh', 'Cdh3', 'Lhx2', 'Sox9', 'Runx1'))

##Cluster 5 marker
FeaturePlot(harmonyEpiCells, features = c('Gli2', 'Bmpr1b', 'Eda', 'Wnt3a', 'Fgfr1', 'Sox6'), ncol = 3)

##Other genes related to hair follicle development
FeaturePlot(harmonyEpiCells, features = c('Gli2', 'Sox9'), split.by = 'genotype')
FeaturePlot(harmonyEpiCells, features = c('Dkk2', 'Lef1'), split.by = 'genotype')
FeaturePlot(harmonyEpiCells, features = c('Eda', 'Edar', 'Nfkb1'), split.by = 'genotype')
FeaturePlot(harmonyEpiCells, features = c('Shh', 'Hoxc13', 'Gli1', 'Msx1', 'Msx2', 'Lgr5', 'Lgr6', 'Pdgfra'))
FeaturePlot(harmonyEpiCells, features = c('Wls', 'Lef1', 'Wnt10b', 'Ctnnb1', 'Dkk4'))


##### Multi-violin plots to highlight select genes #####


#### Define a function to do multiple vlnplot in a column
CB_multi_vlnplot <- function(seurat_cells, features, split.by = NULL, cols = c('darkgreen', 'gold')){
  require(patchwork)
  
  len <- length(features)
  plots <- VlnPlot(seurat_cells, features = features, split.by = split.by, alpha = 0, combine = FALSE,
                   cols = cols)
  plots <- lapply(plots, function(p){
    p <- p + 
      scale_y_continuous(breaks = scales::breaks_pretty(n = 2)) +
      labs(y = p$labels$title) +
      theme(axis.text.x = element_blank(),
            plot.title = element_blank(),
            legend.position = 'none',
            axis.title.x = element_blank(),
            axis.title.y = element_text(angle = 0, vjust = 0.5, size = 12, face = 'bold'),
            plot.margin = margin(r = 20))
  })
  plots[[len]] <- plots[[len]] +
    theme(legend.position = 'right',
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  wrap_plots(plots, ncol = 1)
}

### Genes to be highlighted for aberrDiff
Genes_to_plot <- c('Krt14', 'Krt10', 'Flnb', 'Nectin4', 'Krt17', 'Sox9', 'Gria3', 'Tdrd1', 'Cd28')
CB_multi_vlnplot(harmonyEpiCells, features = Genes_to_plot, split.by = 'genotype')

### Genes to be highlighted (Metabolism, ECM, cytoskeleton) (Export PDF dimension: 4.8 x 6.0)
Genes_to_plot <- c('Pfkl', 'Pgk1', 'Pkm', 'Col28a1', 'Col8a1', 'Lamb1', 'Myh10', 'Tpm1', 'Kank3')
CB_multi_vlnplot(harmonyEpiCells, features = Genes_to_plot, split.by = 'genotype')

### Genes to be highlighted (Cell cycle, p53-related)
Genes_to_plot <- c('Polk', 'Plk2', 'Ddit4l', 'Ccnb1', 'Cdkn1a', 'Ccnd1', 'Wnt3', 'Fgfbp3', 'Ngf')
CB_multi_vlnplot(harmonyEpiCells, features = Genes_to_plot, split.by = 'genotype')

### Genes known to be important for epidermal development
Genes_to_plot <- c('Trp63', 'Klf4', 'Irf6', 'Ovol1', 'Id1', 'Myc', 'Foxn1', 'Gli2', 'Chuk')
CB_multi_vlnplot(harmonyEpiCells, features = Genes_to_plot, split.by = 'genotype')

### Genes that are upregulated in most if not all populations
FeaturePlot(harmonyEpiCells, features = c('Kyat3', 'Tcfl5', 'Plod2'), split.by = 'genotype',
            order = TRUE, by.col = FALSE,
            pt.size = 4, raster = TRUE, raster.dpi = c(300, 300))

### Multi-violin plots of inflammation related genes
Setdb1_cKO <- readRDS('/Users/kb744/Tumbar Lab Local/scRNAseq_Setdb1_cKO_Ana/Seurat_Analysis/Harmony_int_EpithelialCells.RData')

Genes_to_plot <- c('Aim2', 'Ifih1', 'Ddx58', 'Dtx3l', 'Parp9', 'Isg15', 'Creb3', 'Eif2ak2')

p1 <- CB_multi_vlnplot(Setdb1_cKO, features = Genes_to_plot, split.by = 'genotype')
p2 <- CB_multi_vlnplot(harmonyEpiCells, features = Genes_to_plot, split.by = 'genotype')
p1 | p2


### Feature plot for genes highlighting mixed identity
Genes_to_plot <- c('Krt14', 'Krt10', 'Flnb', 'Nectin4', 'Krt17', 'Sox9', 'Gria3', 'Tdrd1', 'Cd28')
FeaturePlot(harmonyEpiCells, features = c('Krt14', 'Nectin4', 'Sox9', 'Gria3'), split.by = 'genotype',
            order = TRUE)


##### Try to make custom multi-violin plots ######

### Define the order of clusters to show for violin plots
## Plot Aberrant populations at last
harmonyEpiCells$Identity <- factor(Idents(harmonyEpiCells),
                                   levels = c('Basal_1', 'Basal_2', 'Basal_3',
                                              'Diff_1', 'Diff_2', 'Diff_3',
                                              'HairFollicle', 'Merkel',
                                              'AberrBasal', "AberrDiff"))
Idents(harmonyEpiCells) <- 'Identity'

## First define a function to make a single plot
CB_custom_vlnplot <- function(seurat_objt, feature = feature, split.by = split.by, color = color) {
  require(ggplot2)
  
  data <- seurat_objt@assays$RNA$scale.data[feature,]
  meta <- seurat_objt@meta.data
  stopifnot(identical(names(data),rownames(meta)))
  data <- cbind(data, meta)
  
  ggplot(data, aes(x = Identity, y = data, fill = .data[[split.by]])) +
    geom_violin(width = 0.7, scale = 'width', position = position_dodge(width = 0.8)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 2)) +
    scale_fill_manual(values = color) +
    labs(x = NULL, y = feature)
}

##Test
CB_custom_vlnplot(harmonyEpiCells, feature = 'Krt14', split.by = 'genotype', color = c('darkgreen', 'gold'))

## Then define function to make combined multiple violin plots
CB_custom_multi_vlnplot <- function(seurat_cells, features, split.by, cols = c('darkgreen', 'gold')){
  require(patchwork)
  
  len = length(features)
  
  plots <- lapply(features, function(g){
    CB_custom_vlnplot(seurat_cells, feature = g, split.by = split.by, color = cols) +
      theme(axis.text.x = element_blank(),
            axis.title.y = element_text(angle = 0, vjust = 0.5, size = 12, face = 'bold'),
            plot.margin = margin(r = 20),
            legend.position = 'none')
  })
  
  plots[[len]] <- plots[[len]] +
    theme(legend.position = 'right',
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  wrap_plots(plots, ncol = 1)
}

### Genes to be highlighted for aberrDiff (Export PDF dimension: 4 x 7.5)
Genes_to_plot <- c('Krt14', 'Krt10', 'Flnb', 'Nectin4', 'Krt17', 'Sox9', 'Gria3', 'Tdrd1', 'Cd28')
CB_custom_multi_vlnplot(harmonyEpiCells, features = Genes_to_plot, split.by = 'genotype')

### Genes to be highlighted (Metabolism, ECM, cytoskeleton)
Genes_to_plot <- c('Pfkl', 'Pgk1', 'Pkm', 'Col28a1', 'Col8a1', 'Lamb1', 'Myh10', 'Tpm1', 'Kank3')
CB_custom_multi_vlnplot(harmonyEpiCells, features = Genes_to_plot, split.by = 'genotype')

### Genes to be highlighted (Cell cycle, p53-related)
Genes_to_plot <- c('Polk', 'Plk2', 'Ddit4l', 'Ccnb1', 'Cdkn1a', 'Ccnd1', 'Wnt3', 'Fgfbp3', 'Ngf')
CB_custom_multi_vlnplot(harmonyEpiCells, features = Genes_to_plot, split.by = 'genotype')

### Genes known to be important for epidermal development
Genes_to_plot <- c('Trp63', 'Klf4', 'Irf6', 'Ovol1', 'Id1', 'Myc', 'Foxn1', 'Gli2', 'Chuk')
CB_custom_multi_vlnplot(harmonyEpiCells, features = Genes_to_plot, split.by = 'genotype')

### Multi-violin plots of inflammation related genes
Setdb1_cKO <- readRDS('/Users/kb744/Tumbar Lab Local/scRNAseq_Setdb1_cKO_Ana/Seurat_Analysis/Harmony_int_EpithelialCells.RData')
Setdb1_cKO$Identity <- Idents(Setdb1_cKO)

Genes_to_plot <- c('Aim2', 'Ifih1', 'Ddx58', 'Dtx3l', 'Parp9', 'Isg15', 'Creb3', 'Eif2ak2')

p1 <- CB_custom_multi_vlnplot(Setdb1_cKO, features = Genes_to_plot, split.by = 'genotype')
p2 <- CB_custom_multi_vlnplot(harmonyEpiCells, features = Genes_to_plot, split.by = 'genotype')
p1 | p2


### Feature plot for genes highlighting mixed identity
Genes_to_plot <- c('Krt14', 'Krt10', 'Flnb', 'Nectin4', 'Krt17', 'Sox9', 'Gria3', 'Tdrd1', 'Cd28')
FeaturePlot(harmonyEpiCells, features = c('Krt14', 'Nectin4', 'Sox9', 'Gria3'), split.by = 'genotype',
            order = TRUE)



##### Other check-up #####
### Function to check if a gene is dysregulated in any of the scRNA-seq clusters
check_DEs <- function(gene, DE_list){
  res <- lapply(DE_list, function(lt){
    gene %in% lt
  })
  return(unlist(res))
}
check_DEs('Trp63', scRNA_DE)
check_DEs('Klf4', scRNA_DE)
check_DEs('Irf6', scRNA_DE)
check_DEs('Ovol1', scRNA_DE)
check_DEs('Id1', scRNA_DE)
check_DEs('Myc', scRNA_DE)
check_DEs('Foxn1', scRNA_DE)
check_DEs('Gli2', scRNA_DE)
check_DEs('Chuk', scRNA_DE)

Idents(harmonyEpiCells) <- "geno_cluster"
Cluster_DE_FCs <- list()
clusters <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
for (name in clusters) {
  Ident_1 <- paste('TKO', name, sep = '_')
  Ident_2 <- paste('Ctrl', name, sep = '_')
  cat('Finding differentially expressed genes for cluster:', name, '\n')
  Cluster_DE_FCs[[name]] <- FindMarkers(harmonyEpiCells, ident.1 = Ident_1, ident.2 = Ident_2, logfc.threshold = 0.3, min.pct = 0.25)
}

Idents(harmonyEpiCells) <- 'Identity'
rm(clusters, name, Ident_1, Ident_2)
