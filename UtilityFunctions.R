library(AnnotationHub)
library(GenomicFeatures)
library(biomaRt)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)

##Getting annotation from ensembldb
ah = AnnotationHub()
dbs <- query(ah, pattern = c("EnsDb", "Mus Musculus", 110))
Ens_mm39_110 <- dbs[['AH113713']]

#saveRDS(Ens_mm39_110, file = 'DataFiles/Enembldb_mm39_110.rds')
rm(ah, dbs)

#####Make Txdb object from Ensembl
txdb_ensembl <- makeTxDbFromEnsembl(organism = 'Mus musculus', release = 110)

###Setup for BiomaRt
enbl <- useEnsembl(biomart = 'genes', dataset = 'mmusculus_gene_ensembl')


###Set up ggplot themes
library(tidyverse)

theme_set(theme_minimal() + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  axis.line = element_line(color = 'black', linewidth = 1),
                  axis.ticks = element_line(color = 'black'),
                  text = element_text(size = 12, face = 'bold', family = 'Arial'),
                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                  axis.title = element_text(size = 14)
            )
)

genoColor <- c('darkgreen', 'gold')
regColor <- c('up' = "#F8766D", 'ns' = "grey", 'down' = "#00BFC4")

###Export DEseq2 results
write_DEseq_result <- function(res, filename) {
  require(openxlsx)
  
  wb <- createWorkbook()
  addWorksheet(wb, 'All')
  addWorksheet(wb, 'DEgenes')
  
  res$padj <- ifelse(is.na(res$padj), 1, res$padj)
  writeData(wb, 'All', as.data.frame(res), rowNames = TRUE)
  
  DEgenes <- res[res$padj<0.05 & abs(res$log2FoldChange)>1, ]
  writeData(wb, 'DEgenes', as.data.frame(DEgenes), rowNames = TRUE)
  
  saveWorkbook(wb, filename, overwrite = TRUE)
}

###Add annotation info
addGeneInfo <- function(filename, ensembl, ID_type){
  require(AnnotationDbi)
  require(openxlsx)
  sheet_names <- getSheetNames(filename)
  new_filename <- paste0(tools::file_path_sans_ext(filename), '_annotated.xlsx')
  
  wb <- createWorkbook()
  for (name in sheet_names) {
    cat('Processing sheet: ', name, '\n')
    
    df <- read.xlsx(filename, sheet = name)
    colnames(df)[1] <- 'GeneID'
    new_info <- AnnotationDbi::select(ensembl, keys = df$GeneID, keytype = ID_type, 
                                      columns = c('GENEID', 'GENEBIOTYPE', 'DESCRIPTION'))
    colnames(new_info)[1] <- 'GeneID'
    df <- merge(df, new_info, by = 'GeneID', all.x = TRUE)
    
    addWorksheet(wb, name)
    writeData(wb, name, df)
  }
  saveWorkbook(wb, file = new_filename, overwrite = TRUE)
  cat('Done!\n')
}

###Export bed files
CB_export_bed <- function(gr, fileName, identifier = 'Symbol') {
  df <- as.data.frame(gr)
  df$start <- df$start - 1
  df$Score <- 0
  #Customize which fields to include in the output here
  bed_df <- df[, c('seqnames', 'start', 'end', identifier, 'Score', 'strand')]
  write.table(bed_df, file = fileName, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
