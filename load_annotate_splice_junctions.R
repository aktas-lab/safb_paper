library(data.table)
library(GenomicFeatures)

source('./utils.R')
source('./config.R')


sj.files <- dir('data/splice_junctions/', full.names = TRUE)

sample.names <- sub("\\..*", "", basename(sj.files))
samplesheet <- data.table(sample = sample.names,
                          sj.tab = sj.files)

sj.dt <- loadSpliceJunctions(samplesheet, hsa.txdb)
sj.annot.dt <-
  annotateSpliceJunctions(sj.dt       = sj.dt,
                          reps.gr     = hsa.rep.gr,
                          txdb        = hsa.txdb,
                          ensg2symbol = hsa.ensg2symbol,
                          ensg2type   = hsa.ensg2type,
                          peaks.gr    = hsa.peaks.reduced.gr)
system('mkdir output')
write.table(sj.annot.dt, 'output/safb_annotated_junctions.txt', quote = F, row.names = F, sep = '\t')
