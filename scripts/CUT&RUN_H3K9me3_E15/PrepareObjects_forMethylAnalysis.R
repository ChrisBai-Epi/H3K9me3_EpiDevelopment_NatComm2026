library(GenomicFeatures)

#Get and define gene region, promoter region and flanking regions
txdb_ensembl <- makeTxDbFromEnsembl(organism = 'Mus musculus', release = 110)
all_genes <- genes(txdb_ensembl)
info <- AnnotationDbi::select(Ens_mm39_110, keys = unlist(all_genes$gene_id), keytype = 'GENEID', 
               columns = c('GENEBIOTYPE', 'GENEID', 'GENENAME', 'SYMBOL', 'DESCRIPTION'))
all_genes$Biotype <- info$GENEBIOTYPE
all_genes$Symbol <- info$SYMBOL
rm(info)

seqlevelsStyle(all_genes) <- 'UCSC'
chrs <- c(paste0('chr', 1:19), 'chrX', 'chrY', 'chrM')
all_genes <- keepSeqlevels(all_genes, chrs, pruning.mode = 'tidy')
all_genes <- all_genes[order(seqnames(all_genes), start(all_genes))]
names(all_genes) <- NULL

all_promoters_10k <- promoters(all_genes, upstream = 10000, downstream = 10000, use.names = TRUE)
