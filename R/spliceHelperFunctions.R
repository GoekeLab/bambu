
#'@title MYGAPS
#'@param x
#'@param start
#'@param end
#'@importFrom GenomicRanges
myGaps <- function(x, start=NA, end=NA)
{
  # License note: This function is adopted from the GenomicAlignments package (Author: Hervé Pagès, Valerie Obenchain, Martin Morgan)
  # License	Artistic-2.0
  # https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments
  if (!.isNumericOrNAs(start))
    stop("'start' must be an integer vector or NA")
  if (!is.integer(start))
    start <- as.integer(start)
  if (!.isNumericOrNAs(end))
    stop("'end' must be an integer vector or NA")
  if (!is.integer(end))
    end <- as.integer(end)

  ## seqname and strand consistent in list elements
  if (all(elementNROWS(runValue(seqnames(x))) == 1L) &&
      all(elementNROWS(runValue(strand(x))) == 1L)) {
    flat <- unlist(x, use.names=FALSE)
    gaps <- gaps(ranges(x), start, end)
    ### FIXME: this makes this function more of an 'introns' than a .gaps.
    ### FIXME: this breaks when the GRangesList is not ordered by position
    if (!is.null(mcols(x, use.names=FALSE)$query.break)) {
      insert_gaps <- as(ranges(.insertGaps(x)), "CompressedIRangesList")
      gaps <- setdiff(gaps, insert_gaps)
    }

    idx <- elementNROWS(gaps) != 0
    ## FIXME : can't handle lists with empty elements
    ##         'start' and 'end' not quite right here
    firstseg <- start(PartitioningByWidth(x))
    seqnms <- rep(seqnames(flat)[firstseg], elementNROWS(gaps))
    strand <- rep(strand(flat)[firstseg],elementNROWS(gaps))
    gr <- relist(GRanges(seqnms, unlist(gaps, use.names=FALSE), strand), gaps)
    gr
  } else {
    ### FIXME: does not handle query.break column yet
    setdiff(range(x), x)
  }

}

.isNumericOrNAs <- S4Vectors:::isNumericOrNAs






##### UNTIL HERE ##



#'@title RANGEDIST
#'@param query
#'@param subject
#'@param splice
#'@param maxDist
rangesDist <- function(query, subject, splice, maxDist)
{
  qrng <- ranges(query)
  srng <- ranges(subject)
  sprng <- ranges(splice)

  setDiffQ <-   width(GenomicRanges::setdiff(qrng, srng))
  interesectS <- width(GenomicRanges::intersect(srng, sprng))
  uniqueLengthQuery <- sum(setDiffQ)
  uniqueLengthSubject <- sum(interesectS)

  queryElementsOutsideMaxDist <- sum(setDiffQ>=maxDist)
  subjectElementsOutsideMaxDist <-  sum(interesectS>=maxDist)
  compatible <- (queryElementsOutsideMaxDist == 0) & (subjectElementsOutsideMaxDist == 0)
  DataFrame(uniqueLengthQuery,uniqueLengthSubject, compatible, queryElementsOutsideMaxDist, subjectElementsOutsideMaxDist) ##### HERE ####
}

#'@title FINDSPLICEOVERLAPSQUICK
#'@description the following functions are implemented in R, I just included the within option to make them significantly faster+memorey friendly for this purpose (original code copied from https://rdrr.io/bioc/GenomicAlignments/src/R/findSpliceOverlaps-methods.R)
#'@param query
#'@param subject
#'@param ignore.strand
findSpliceOverlapsQuick <- function(query, subject, ignore.strand=FALSE) {
  olap <- findOverlaps(query, subject, ignore.strand=ignore.strand, type='within')
  olapEqual <- findOverlaps(query, subject, ignore.strand=ignore.strand, type='equal')
  if (length(olap) == 0L)
    return(.result(olap))


  query <- query[queryHits(olap)]
  subject <- subject[subjectHits(olap)]
  splice <- myGaps(query)


  compatible <- myCompatibleTranscription(query, subject, splice)
  strandSpecific <- all(strand(query) != "*")
  rm(list=c('query','subject'))
  gc()
  equal <- (!is.na(S4Vectors::match(olap ,olapEqual)))
  rm(olapEqual)
  gc()
  unique <- myOneMatch(compatible, queryHits(olap))

  mcols(olap) <- DataFrame(compatible, equal, unique, strandSpecific)
  olap
}

#'@title MYCOMPATIBLETRANSCRIPTION
#'@param splice
#'@param query
#'@param subject
myCompatibleTranscription <- function(query, subject, splice)
{
  qrng <- ranges(query)
  srng <- ranges(subject)
  sprng <- ranges(splice)

  bnds <- elementNROWS(GenomicRanges::setdiff(qrng, srng)) == 0L
  splc <- elementNROWS(GenomicRanges::intersect(srng, sprng)) == 0L
  return(bnds & splc)
}

#'@title MYONEMATCH
#'@param x
#'@param idx
myOneMatch <- function(x, idx)
{
  # License note: This function is adopted from the GenomicAlignments package (Author: Hervé Pagès, Valerie Obenchain, Martin Morgan)
  # https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments
  xcnt <- rowsum(as.integer(x), idx)[,1]
  oneMatch <- rep((xcnt == 1L), table(idx))
  unname(x & oneMatch)
}



#'@title .EXTRACT_UNORIENTED_INTRON_MOTIF
#'@param genome
#'@param junctions
.extract_unoriented_intron_motif <- function(genome, junctions)
{
  mcols(junctions) <- NULL
  junctions_len <- length(junctions)
  Ldinucl_gr <- Rdinucl_gr <- junctions
  end(Ldinucl_gr) <- start(Ldinucl_gr) + 1L
  start(Rdinucl_gr) <- end(Rdinucl_gr) - 1L
  all_dinucl <- getSeq(genome, c(Ldinucl_gr, Rdinucl_gr))
  Rdinucl <- tail(all_dinucl, n=junctions_len)
  xscat(Ldinucl, "-", Rdinucl)
}

#'@title SPLICESTRAND
#'@param motif
spliceStrand <- function(motif){
  NATURAL_INTRON_MOTIFS_RC <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(GenomicAlignments::NATURAL_INTRON_MOTIFS)))

  motifStrand <- ifelse(motif %in% GenomicAlignments::NATURAL_INTRON_MOTIFS,'+','*')
  motifStrand[motif %in% NATURAL_INTRON_MOTIFS_RC] <- '-'
  return(motifStrand)
}


#'@title TRIMFIRSTLASTEXONS
#'@description function returns exon granges list where the first and last exons are trimmed to width 2, all introns are preserved, requires exon_rank and exon_endRank
#'@param grangesListWithExonRanks
trimFirstLastExons <- function(grangesListWithExonRanks) {
  unlistedGrangesList <- unlist(grangesListWithExonRanks, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(grangesListWithExonRanks)), names=NULL)
  startExonsSet <- (which((unlistedGrangesList$exon_rank==1 & as.character(strand(unlistedGrangesList))!='-')|(unlistedGrangesList$exon_endRank==1 & as.character(strand(unlistedGrangesList))=='-')))
  endExonsSet <- (which((unlistedGrangesList$exon_rank==1 & as.character(strand(unlistedGrangesList))=='-')|(unlistedGrangesList$exon_endRank==1 & as.character(strand(unlistedGrangesList))!='-')))
  start(unlistedGrangesList[startExonsSet]) <-  end(unlistedGrangesList[startExonsSet])-1
  end(unlistedGrangesList[endExonsSet]) <-  start(unlistedGrangesList[endExonsSet])+1
  return(relist(unlistedGrangesList, partitioning))
}

#'@title MYPERFORMANCE
#'@param labels
#'@param scores
#'@param descreasing
myPerformance <- function(labels, scores, decreasing = TRUE){
  labels <- labels[order(scores, decreasing = decreasing)]
  results <- list()
  results[['TPR']] <- cumsum(labels)/sum(labels)  # TP/(TP+FP); True Positive Rate;Sensitivity; recall
  results[['FPR']] <- cumsum(!labels)/sum(!labels)  # FP/(FP+TN); False Positive Rate;1-Specificity
  results[['precision']] <- cumsum(labels)/(1:length(labels))  # TP/(TP+FP); positive predictive value;precision
  results[['AUC']] <- sum(results[['TPR']][!duplicated( results[['FPR']],fromLast=TRUE)]/sum(!duplicated( results[['FPR']],fromLast=TRUE)))
  return(results)
}


#'@title PLOTGRANGESLISTBYNAMES
#'@param grangesList
plotGRangesListByNames<-function(grangesList)
{

  plotTrack = unlist(grangesList)
  plotTrack$grouping <- names(plotTrack)
  plotTrack=AnnotationTrack(plotTrack,group=plotTrack$grouping, id=plotTrack$grouping,fill=2,col.line=1,fontcolor.group=1,fontcolor.item=1)
  plotTracks(list(plotTrack),showId=TRUE)
}

#' @describeIn plotGRangesListByNames
makeTrackFromGrangesList <- function(grangesList)
{

  plotTrack = unlist(grangesList)
  plotTrack$grouping <- names(plotTrack)
  plotTrack=AnnotationTrack(plotTrack,group=plotTrack$grouping, id=plotTrack$grouping,fill=2,col.line=1,fontcolor.group=1,fontcolor.item=1)
  return(plotTrack)
}


#' @title GRANGESLISTTOBED
#' @param x
grangesListToBed<-function(x) # note: 0 based coordinates
{
  # useful: generate UCSC custom track with the following steps


  xUnlist <- unlist(range(x))

  bedData = data.frame(as.character(seqnames(xUnlist)),start(xUnlist)-1,end(xUnlist),names(xUnlist),rep(1000,length(xUnlist)),as.character(strand(xUnlist)),start(xUnlist)-1,end(xUnlist),rep(0,length(xUnlist)),elementNROWS(x),unlist(lapply(width(x),paste,collapse=',')),unlist(lapply((start(x)-min(start(x))),paste,collapse=','))) #, stringsAsFactors=FALSE
  colnames(bedData) <-  c('chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts')
  return(bedData)
}
