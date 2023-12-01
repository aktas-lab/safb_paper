hsa.gtf.gr <- readRDS('./data/annotation/gencode.human.v29.rds')

hsa.ensg2type <- unique(as.data.table(hsa.gtf.gr)[type == 'gene'][,.(gene_id, gene_type)][!is.na(gene_id)])
hsa.ensg2symbol <- unique(as.data.table(mcols(hsa.gtf.gr)[, c('gene_id', 'gene_name')]))
setnames(hsa.ensg2symbol, c('gene_id', 'symbol'))

hsa.txdb <- makeTxDbFromGRanges(hsa.gtf.gr)
hsa.genes.gr <-
  keepStandardChromosomes(genes(hsa.txdb),
                          pruning.mode = 'coarse')
hsa.rep.gr <-
  keepStandardChromosomes(
    readRDS("./data/annotation/hg38.repeats.classes.flipped.rds"),
    pruning.mode = 'coarse')

hsa.peaks.reduced.gr <-
  readRDS('./data/annotation/human_SAFB_SAFB2_SLTM_peaks_reduced_gr.rds')