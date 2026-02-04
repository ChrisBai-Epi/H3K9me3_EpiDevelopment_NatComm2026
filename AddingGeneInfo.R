library(openxlsx)
library(biomaRt)

filename <- "ExcelFiles/Conserved_markers_v3.xlsx"
enbl <- useEnsembl(biomart = 'genes', dataset = 'mmusculus_gene_ensembl')

addGeneInfo <- function(filename, ensembl){
  require(biomaRt)
  require(openxlsx)
  sheet_names <- getSheetNames(filename)
  new_filename <- paste0(tools::file_path_sans_ext(filename), '_annotated.xlsx')
  
  wb <- createWorkbook()
  for (name in sheet_names) {
    cat('Processing sheet: ', name, '\n')
    
    df <- read.xlsx(filename, sheet = name)
    colnames(df)[1] <- 'GeneID'
    new_info <- getBM(mart = enbl, attributes = c('external_gene_name', 'wikigene_description', 'gene_biotype'), 
                      filters = 'external_gene_name', df$GeneID)
    colnames(new_info)[1] <- 'GeneID'
    new_info <- new_info[!duplicated(new_info$GeneID),]
    df <- merge(df, new_info, by = 'GeneID', all.x = TRUE)
    
    addWorksheet(wb, name)
    writeData(wb, name, df)
  }
  saveWorkbook(wb, file = new_filename, overwrite = TRUE)
  cat('Done!\n')
}

addGeneInfo("ExcelFiles/Conserved_markers_v3.xlsx", enbl)
addGeneInfo("ExcelFiles/DE_genes_v3.xlsx", enbl)
addGeneInfo("ExcelFiles/TKO_all_vs_Ctrl_all.xlsx", enbl)
addGeneInfo("ExcelFiles/TKO_basal_vs_Ctrl_basal.xlsx", enbl)
addGeneInfo("ExcelFiles/TKO_diff_vs_Ctrl_diff.xlsx", enbl)
