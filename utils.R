readSJ <- function(filename) {
  # reads STAR SJ.out.tab
  
  sj.dt <- fread(filename)
  setnames(sj.dt, c("chr", "start", "end", "strand",
                    "motif", "known",
                    "uniq.mappers", "multimappers", "max.overhang"))
  sj.dt[, pos := paste0(chr, ":", start, "-", end)]
  
  sj.dt
  
}

loadSpliceJunctions <- function(star.files, txdb) {
  # loads and annotates STAR splice junctions
  
  sj.ls <- parallel::mclapply(star.files$sj.tab, readSJ, mc.cores = 6)
  names(sj.ls) <- star.files$sample
  
  sj.summ.dt <- summarizeList(sj.ls, "pos", "uniq.mappers")
  sj.summ.dt[is.na(sj.summ.dt)] <- 0
  sj.summ.dt <- sj.summ.dt[rowSums(sj.summ.dt[, grep("rep", names(sj.summ.dt)), with = F]) > 0]
  
  introns.gr <- unique(unlist(intronsByTranscript(txdb)))
  introns.gr <- keepStandardChromosomes(introns.gr, pruning.mode = 'coarse')
  genes.gr  <- genes(txdb)
  genes.gr <- keepStandardChromosomes(genes.gr, pruning.mode = 'coarse')
  
  sj.annot.dt <- rbindlist(sj.ls)[,.(chr, start, end, motif2 = motif, strand2 = strand, known)]
  sj.annot.dt <- unique(sj.annot.dt)
  sj.annot.dt <- sj.annot.dt[strand2 != 0]
  sj.annot.dt <- sj.annot.dt[known == 1]
  sj.annot.dt[, strand := ifelse(strand2 == 1, "+", "-")]
  
  sj.annot.dt[, motif := '']
  sj.annot.dt[motif2 == 0, motif := 'noncanonical']
  sj.annot.dt[motif2 == 1, motif := 'GT/AG']
  sj.annot.dt[motif2 == 2, motif := 'CT/AC']
  sj.annot.dt[motif2 == 3, motif := 'GC/AG']
  sj.annot.dt[motif2 == 4, motif := 'CT/GC']
  sj.annot.dt[motif2 == 5, motif := 'AT/AC']
  sj.annot.dt[motif2 == 6, motif := 'GT/AT']
  
  # drop splices from non-standard chromosomes
  sj.annot.dt <- sj.annot.dt[chr %in% seqlevels(genes.gr)]
  
  sj.annot.gr <- dt2gr(sj.annot.dt[,.(chr, start, end, strand)])
  
  # these should be annotated introns
  olaps.equal <- findOverlaps(sj.annot.gr, introns.gr, type = "equal")
  
  # annotated left coordinate
  olaps.left  <- findOverlaps(sj.annot.gr, introns.gr, type = "start")
  
  # annotated right coordinate
  olaps.right <- findOverlaps(sj.annot.gr, introns.gr, type = "end")
  
  # within a gene?
  olaps.within.a.gene <- findOverlaps(sj.annot.gr, genes.gr, type = "within")
  
  # start within a gene
  sj.start.annot.gr <- resize(sj.annot.gr, 1, fix = "start")
  sj.end.annot.gr   <- resize(sj.annot.gr, 1, fix = "end")
  olaps.start.within.a.gene <- findOverlaps(sj.start.annot.gr, genes.gr, type = "within")
  olaps.end.within.a.gene <- findOverlaps(sj.end.annot.gr, genes.gr, type = "within")
  
  sj.annot.dt[, idx := 1:nrow(sj.annot.dt)]
  sj.annot.dt[, junct.known := "none"]
  sj.annot.dt[strand == "+" & idx %in% queryHits(olaps.left),  junct.known := "don"]
  sj.annot.dt[strand == "+" & idx %in% queryHits(olaps.right), junct.known := "acc"]
  sj.annot.dt[strand == "-" & idx %in% queryHits(olaps.left),  junct.known := "acc"]
  sj.annot.dt[strand == "-" & idx %in% queryHits(olaps.right), junct.known := "don"]
  
  sj.annot.dt[idx %in% queryHits(olaps.left) & idx %in% queryHits(olaps.right), junct.known := "both"]
  sj.annot.dt[idx %in% queryHits(olaps.equal), junct.known := "annot"]
  
  sj.annot.dt[, ingene := idx %in% queryHits(olaps.within.a.gene)]
  
  sj.annot.dt[, flankInGene := idx %in% queryHits(olaps.start.within.a.gene) | idx %in% queryHits(olaps.end.within.a.gene)]
  
  sj.annot.dt[, pos := paste0(chr, ":", start, "-", end)]
  
  
  sj.out.dt <- merge(sj.annot.dt[,.(pos, strand, motif, junct.known, ingene)],
                     sj.summ.dt,
                     by = "pos")
  
  return(sj.out.dt)
}

annotateSpliceJunctions <- function(sj.dt, reps.gr, txdb, ensg2symbol, ensg2type, peaks.gr) {
  # sj.dt - STAR SJ.out.tab preprocessed with loadSpliceJunctions()
  # reps.gr - repeat annotation with classes, families,
  #           and antisense elements
  # txdb - transcript database
  # ensg2symbol - ensembl gene id xrefs to gene symbols
  # txdb.pseudo - pseudogene txdb. if !hasArg, skip

  # prep return data.table
  sj.out.dt <- data.table::copy(sj.dt)
  
  # split position into chr/start/end
  # prep GR with splice donor and acceptor position
  DT <- data.table()
  DT[, c('chr', 'start', 'end') := tstrsplit(sj.out.dt$pos, ":|-")]
  DT[, start := as.integer(start)]
  DT[, end := as.integer(end)]
  DT[, strand := sj.out.dt$strand]
  
  
  sj.gr <- dt2gr(DT)
  sj.start.gr <- resize(sj.gr, 1, fix = 'start')
  sj.end.gr   <- resize(sj.gr, 1, fix = 'end')
  
  # split repeats by family
  reps.family.grl <- split(reps.gr, reps.gr$family)
  
  # are donor and acceptor sitting in a repeat? which one?
  sj.out.dt[, don.rep := annotateRanges(sj.start.gr,
                                                reps.family.grl,
                                                type = 'all',
                                                null.fact = '')]
  
  sj.out.dt[, acc.rep := annotateRanges(sj.end.gr,
                                                reps.family.grl,
                                                type = 'all',
                                                null.fact = '')]
  
  # is the thing getting spliced to a gene?
  # some more annotation prepping
  ensg2enst <-
    unique(
      data.table(
        select(txdb,
               keys = keys(txdb),
               columns = c('GENEID', 'TXNAME'),
               keytype = 'GENEID')))
  setnames(ensg2enst, c('gene_id', 'tx_id'))
  enst2symbol <- merge(ensg2enst, ensg2symbol, 'gene_id', all.x = TRUE)
  
  exbt.grl <- exonsBy(txdb, by = 'tx', use.names = TRUE)
  exons.gr <- unique(unlist(exbt.grl))
  exons.gr <- keepStandardChromosomes(exons.gr, pruning.mode = 'coarse')
  exons.gr$symbol <- enst2symbol[match(names(exons.gr), enst2symbol$tx_id)]$symbol
  
  don.in.exon <- findOverlaps(sj.start.gr, exons.gr, type = 'within')
  acc.in.exon <- findOverlaps(sj.end.gr, exons.gr, type = 'within')
  
  sj.out.dt[queryHits(don.in.exon), don.exon := exons.gr[subjectHits(don.in.exon)]$symbol]
  sj.out.dt[queryHits(acc.in.exon), acc.exon := exons.gr[subjectHits(acc.in.exon)]$symbol]
  
  # is one or both of the splice sites in peaks?
  if (hasArg(peaks.gr)) {

    # labels those inside peaks
    don.peak <- overlapsAny(sj.start.gr,
                            peaks.gr,
                            type = 'within',
                            ignore.strand = TRUE)
    acc.peak <- overlapsAny(sj.end.gr,
                            peaks.gr,
                            type = 'within',
                            ignore.strand = TRUE)
    
    mrg.peak <- rep('none', length(don.peak))
    mrg.peak[don.peak == TRUE  & acc.peak == TRUE]  <- 'both'
    mrg.peak[don.peak == TRUE  & acc.peak == FALSE] <- 'don'
    mrg.peak[don.peak == FALSE & acc.peak == TRUE]  <- 'acc'
    
    sj.out.dt[, ss.peak := mrg.peak]
  }
  
  # annotate genes
  # prep annotation
  genes.gr  <- genes(txdb)
  genes.gr <- keepStandardChromosomes(genes.gr, pruning.mode = 'coarse')
  
  # to avoid one-to-many relationships, split genes by protein-coding and other
  # sort both sets by length descendingly, and concatenate
  # in olaps, drop all duplicated queryHits
  # this way, mapping to protein_coding is prioritized over other genes
  # and longer genes are prioritized over shorter genes
  genes.pc.gr <- genes.gr[names(genes.gr) %in% ensg2type[gene_type == 'protein_coding']$gene_id]
  genes.pc.gr <- genes.pc.gr[order(-width(genes.pc.gr))]
  
  genes.nonpc.gr <- genes.gr[!(names(genes.gr) %in% ensg2type[gene_type == 'protein_coding']$gene_id)]
  genes.nonpc.gr <- genes.nonpc.gr[order(-width(genes.nonpc.gr))]

  genes.sorted.gr <- c(genes.pc.gr, genes.nonpc.gr)
  
  golaps <- findOverlaps(sj.gr, genes.sorted.gr, type = 'within')
  golaps <- golaps[!duplicated(queryHits(golaps))]
  
  sj.out.dt[queryHits(golaps), gene_id := genes.sorted.gr[subjectHits(golaps)]$gene_id]
  sj.out.dt <- merge(sj.out.dt, ensg2symbol, by = 'gene_id', all.x = TRUE)
  
  return(sj.out.dt)
}

annotateRanges <- function(r1, l,
                           ignore.strand = FALSE,
                           type = "precedence",
                           null.fact = "None",
                           collapse.char = ":") {
    if (!class(r1) == "GRanges")
        stop("Ranges to be annotated need to be GRanges")
    if (!all(sapply(l, class) == "GRanges"))
        stop("Annotating ranges need to be GRanges")
    if (!type %in% c("precedence", "all"))
        stop("type may only be precedence and all")
    if (class(l) != "GRangesList")
        l <- GRangesList(lapply(l, function(x) {
            values(x) <- NULL
            x
        }))
    a <- suppressWarnings(data.table(as.matrix(findOverlaps(r1,
        l, ignore.strand = ignore.strand))))
    a$id <- names(l)[a$subjectHits]
    a$precedence <- match(a$id, names(l))
    a <- a[order(a$precedence)]
    if (type == "precedence") {
        a <- a[!duplicated(a$queryHits)]
    }
    if (type == "all") {
        a <- a[, list(id = paste(unique(id), collapse = collapse.char)),
            by = "queryHits"]
    }
    annot <- rep(null.fact, length(r1))
    annot[a$queryHits] <- a$id
    return(annot)
}

summarizeList <- function (LS, id, value) {
    tmp.dt <- lapply(LS, function(x) x[, c(id, value), with = F])
    lapply(names(tmp.dt), function(x) setnames(tmp.dt[[x]], 2,
        paste(x, "value", sep = ".")))
    DT <- Reduce(function(x, y) merge(x, y, by = id, all = TRUE),
        tmp.dt)
    setnames(DT, sub(".value", "", names(DT)))
    DT
}

dt2gr <- function(DT) {
  gr <- GRanges(seqnames = DT[[1]],
                ranges   = IRanges(start = DT[[2]],
                                   end   = DT[[3]]),
                strand   = DT[[4]])
  # add metadata if exists
  if (ncol(DT) > 4) {
    values(gr) <- DT[,5:ncol(DT)]
  }

  return(gr)
}