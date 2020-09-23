# Note: several of the functions in this file are adopted from the GenomicAlignments package (Author: Hervé Pagès, Valerie Obenchain, Martin Morgan)
# License	Artistic-2.0
# https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments


#' This function calcualtes compatible splice overlaps allowing for a distance threshold, and returns distance in bp between query and subject. Can be used to assign more transcripts to annotations and reads to transcripts.
#' @noRd
findSpliceOverlapsByDist <-function(query, subject, ignore.strand=FALSE, maxDist = 5, type='within', firstLastSeparate = TRUE, dropRangesByMinLength=FALSE, cutStartEnd = TRUE) {

  #  with this option the first and last exons are stored and the distance for each between query and subject hits is returned
  if(firstLastSeparate) {
    queryStart <- selectStartExonsFromGrangesList(query, exonNumber = 1)
    queryEnd <- selectEndExonsFromGrangesList(query, exonNumber = 1)
    subjectStart <- selectStartExonsFromGrangesList(subject, exonNumber = 1)
    subjectEnd <- selectEndExonsFromGrangesList(subject, exonNumber = 1)
    subjectFull <- subject
  }


  #this list is used to find candidate matches. allows ranges to be dropped from query (to find matches not strictly within). more matches>slower
  if(dropRangesByMinLength) {
    queryForOverlap <- dropGrangesListElementsByWidth(query, minWidth=maxDist, cutStartEnd=cutStartEnd)
  } else if(cutStartEnd) {
    queryForOverlap <- cutStartEndFromGrangesList(query)
  } else {
    queryForOverlap <- query
  }


  # first and last exons are reduced to 2bp, internal exon and splice matches are emphasised
  query <- cutStartEndFromGrangesList(query)
  #  subject <- cutStartEndFromGrangesList(subject)

  #  queryExonsRemoved <- dropGrangesListElementsByWidth(query, minWidth=maxDist)
  subjectExtend <- extendGrangesListElements(subject, by=maxDist)

  olap <- findOverlaps(queryForOverlap, subjectExtend, ignore.strand=ignore.strand, type=type)
  olapEqual <- findOverlaps(query, cutStartEndFromGrangesList(subject), ignore.strand=ignore.strand, type='equal')

  query <- query[queryHits(olap)]
  subject <- subject[subjectHits(olap)]
  splice <- myGaps(query)


  compatible <- rangesDist(query, subject, splice, maxDist)
  equal <- (!is.na(S4Vectors::match(olap, olapEqual)))
  unique <- myOneMatch(compatible$compatible, queryHits(olap))
  strandSpecific <- all(strand(query) != "*")
  strandedMatch <- ((all(strand(query)=='-') & all(strand(subject)=='-'))| (all(strand(query)=='+') & all(strand(subject)=='+')))
  mcols(olap) <- DataFrame(compatible, equal, unique, strandSpecific, strandedMatch)

  if(firstLastSeparate) { ##NOTE: Check if there is an error with the start sequence ##
    if(length(olap)>0) {
      queryStart <- ranges(queryStart[queryHits(olap)])
      subjectStart <- ranges(subjectStart[subjectHits(olap)])
      queryEnd <- ranges(queryEnd[queryHits(olap)])
      subjectEnd <- ranges(subjectEnd[subjectHits(olap)])
      subjectFull <- ranges(subjectFull[subjectHits(olap)])
      #queryFirstJunction <- ranges(selectStartExonsFromGrangesList(splice, exonNumber = 1))
      #queryLastJunction <- ranges(selectEndExonsFromGrangesList(splice, exonNumber = 1))

      startMatch <- poverlaps(unlist(queryStart), unlist(subjectStart))
      endMatch <- poverlaps(unlist(queryEnd), unlist(subjectEnd))

      #calculate distance between first and last exon matches:
      # distances corresponds to length of unique sequences in all matching exons
      # distance = width(element in query) + width(all matching elements in subject) - 2x width(intersection(query, subject)
      subjectList <- unlist(subjectFull)
      queryStartList <- rep(unlist(queryStart),elementNROWS(subjectFull))
      myId <- rep(1:length(queryStart),elementNROWS(subjectFull))
      byExonIntersect=pintersect(queryStartList,subjectList, resolve.empty='start.x')
      startDist <- as.integer(tapply(width(subjectList)*(width(byExonIntersect)>0)-2*width(byExonIntersect),myId,sum) + width(unlist(queryStart)))
      uniqueStartLengthQuery <-   sum(width(GenomicRanges::setdiff(queryStart, subjectFull)))
      uniqueStartLengthSubject <-  startDist - uniqueStartLengthQuery



      queryEndList <- rep(unlist(queryEnd),elementNROWS(subjectFull))
      myId <- rep(1:length(queryEnd),elementNROWS(subjectFull))
      byExonIntersect=pintersect(queryEndList,subjectList, resolve.empty='start.x')
      endDist <- as.integer(tapply(width(subjectList)*(width(byExonIntersect)>0)-2*width(byExonIntersect),myId,sum) + width(unlist(queryEnd)))
      uniqueEndLengthQuery <-   sum(width(GenomicRanges::setdiff(queryEnd, subjectFull)))
      uniqueEndLengthSubject <-  endDist - uniqueEndLengthQuery
    } else {
      startMatch <- NULL
      uniqueStartLengthQuery <- NULL
      uniqueStartLengthSubject <- NULL
      endMatch <- NULL
      uniqueEndLengthQuery <- NULL
      uniqueEndLengthSubject <- NULL
    }

    mcols(olap) <- DataFrame( mcols(olap), startMatch, uniqueStartLengthQuery,uniqueStartLengthSubject,endMatch, uniqueEndLengthQuery,uniqueEndLengthSubject)
  }

  olap
}

#' annotate splice overlap by distance
#' @noRd
annotateSpliceOverlapsByDist <-function(query, subject) {
  subjectFullRng <- ranges(subject)
  queryFullRng <- ranges(query)
  strand <- as.character(getStrandFromGrList(query))
  queryStartRng <- selectStartExonFromRangesList(queryFullRng, strand)
  subjectStartRng <- selectStartExonFromRangesList(subjectFullRng, strand)
  queryEndRng <- selectEndExonFromRangesList(queryFullRng, strand)
  subjectEndRng <- selectEndExonFromRangesList(subjectFullRng, strand)
  querySpliceRng <- ranges(myGaps(query))
  querySpliceRng[elementNROWS(querySpliceRng)==0] <- IRanges(start=1,end=1) # add mock intron
  subjectSpliceRng <- ranges(myGaps(subject))
  subjectSpliceRng[elementNROWS(subjectSpliceRng)==0] <- IRanges(start=1,end=1)# add mock intron
  ###############################################################################
  annotatedTable <- tibble(queryId = names(query), subjectId = names(subject), strand = strand)
  # calculate alternative First/last exons and annotate internal start and end first exons
  alterEndTable <- annotateAlterEndExons(queryStartRng, queryEndRng, queryFullRng,
                                         subjectStartRng, subjectEndRng, subjectFullRng, annotatedTable$strand)
  #intron retention subject
  annotatedTable$intronRetention.subject <- annotateIntronRetent(querySpliceRng, subjectFullRng)
  #intron retention query
  annotatedTable$intronRetention.query <- annotateIntronRetent(subjectSpliceRng ,queryFullRng)
  #exon skipping query
  annotatedTable$exonSkipping.query <- annotateExonSkip(querySpliceRng, subjectFullRng,
                                                        subjectStartRng, subjectEndRng)
  #exon skipping subject
  annotatedTable$exonSkipping.subject <- annotateExonSkip(subjectSpliceRng ,queryFullRng,
                                                          queryStartRng, queryEndRng) 
  # exon 3' splice site
  exonSpliceTable <- annotateExonSplice(querySpliceRng, subjectFullRng,
                                        subjectStartRng, subjectEndRng,
                                        annotatedTable$strand)
  annotatedTable <- cbind(annotatedTable, alterEndTable, exonSpliceTable)
  return(annotatedTable)
}

#' extract strand from GRangesList
#' @noRd
getStrandFromGrList <- function(grl) { 
  return(unlist(strand(grl), use.names = F)[cumsum(elementNROWS(grl))]) 
  }

#' start/end ranges pre-processing                                 
expandStartEndRanges <- function(unexpandedRng,targetRng){ 
  processedRng <- rep(unexpandedRng,rep(elementNROWS(targetRng),times=rep(1,length(unexpandedRng))))
  mcols(processedRng)$IdMap <- rep(1:length(unexpandedRng),rep(1,length(unexpandedRng))*elementNROWS(targetRng))
  mcols(processedRng)$matchRng <- unlist(rep(targetRng, times= rep(1,length(unexpandedRng))), use.names=F)
  return (processedRng)
}

#' splice ranges pre-processing
expandSpliceRanges <- function(unexpandedRng,targetRng){ 
  processedRng <- rep(unlist(unexpandedRng, use.names=F),rep(elementNROWS(targetRng),times=elementNROWS(unexpandedRng)))
  mcols(processedRng)$IdMap <- rep(1:length(unexpandedRng),elementNROWS(unexpandedRng)*elementNROWS(targetRng))
  mcols(processedRng)$matchRng <- unlist(rep(targetRng, times= elementNROWS(unexpandedRng)), use.names=F)
  return (processedRng)
}

#' annotate alternative First/last exons and TSS/TES                                  
#' @noRd
annotateAlterEndExons <- function(queryStartRng, queryEndRng, queryFullRng,
                                  subjectStartRng, subjectEndRng, subjectFullRng, strand){
  alterEndTable <-  tibble(start.first.query = start(queryStartRng),
                           end.first.query = end(queryStartRng),
                           start.last.query = start(queryEndRng),
                           end.last.query = end(queryEndRng),
                           start.first.subject = start(subjectStartRng),
                           end.first.subject = end(subjectStartRng),
                           start.last.subject = start(subjectEndRng),
                           end.last.subject = end(subjectEndRng),
                           strand = strand) 
  alterEndTable <- alterEndTable %>% mutate(alternativeFirstExon=!ifelse(strand!='-',
                                                                         start.first.query <= end.first.subject & end.first.query >= start.first.subject,
                                                                         end.first.query >= start.first.subject & start.first.query <= end.first.subject),
                                            alternativeLastExon=!ifelse(strand!='-',
                                                                        end.last.query >= start.last.subject & start.last.query <= end.last.subject,
                                                                        start.last.query <= end.last.subject & end.last.query >= start.last.subject),
                                            alternativeTSS=ifelse(strand!='-',
                                                                  start.first.subject-start.first.query,
                                                                  end.first.query- end.first.subject )*!alternativeFirstExon,
                                            alternativeTES=ifelse(strand!='-', 
                                                                  end.last.query-end.last.subject,
                                                                  start.last.subject-start.last.query)*!alternativeLastExon)
  alterEndTable <- alterEndTable %>% select(alternativeFirstExon, alternativeLastExon, alternativeTSS, alternativeTES)
  #internal start/end query 
  alterEndTable$internalFirstExon.query <- annotateInternalStartEnd(queryStartRng, subjectFullRng, alterEndTable$alternativeFirstExon)
  alterEndTable$internalLastExon.query <- annotateInternalStartEnd(queryEndRng, subjectFullRng, alterEndTable$alternativeLastExon)
  #internal start/end subject
  alterEndTable$internalFirstExon.subject <- annotateInternalStartEnd(subjectStartRng, queryFullRng, alterEndTable$alternativeFirstExon)
  alterEndTable$internalLastExon.subject <- annotateInternalStartEnd(subjectEndRng, queryFullRng, alterEndTable$alternativeLastExon)
  return (alterEndTable)
}

#' annotate internal start and end first exons                                    
#' @noRd
annotateInternalStartEnd <- function(exonRng, fullRng, alternativeFirstLastExon){
  exon.Full.Rng <- expandStartEndRanges(exonRng, fullRng) 
  #internal start/end
  exonIntersect <- pintersect(exon.Full.Rng,mcols(exon.Full.Rng)$matchRng, resolve.empty='start.x')
  internalStartEndVector <- tapply(width(exonIntersect), mcols(exon.Full.Rng)$IdMap, sum)!= 0 & alternativeFirstLastExon
  return(internalStartEndVector)
}

#' annotate intron retention 
#' @noRd
annotateIntronRetent <- function(spliceRng, fullRng){
  splice.FullSplice.Rng <- expandSpliceRanges(spliceRng, fullRng)
  intronRetention <-punion(splice.FullSplice.Rng, mcols(splice.FullSplice.Rng)$matchRng, fill.gap=TRUE) == mcols(splice.FullSplice.Rng)$matchRng
  intronRetentionVector <- tapply(intronRetention, mcols(splice.FullSplice.Rng)$IdMap, sum)
  return (intronRetentionVector)
}

#' annotate exon skiping 
#' @noRd
annotateExonSkip <- function(spliceRng, fullRng, startRng, endRng){
  splice.FullSplice.Rng <- expandSpliceRanges(spliceRng, fullRng)
  start.Splice.Rng <- expandStartEndRanges(startRng, spliceRng)
  end.Splice.Rng <- rep(endRng, elementNROWS(spliceRng))
  exonSkipping <- punion(splice.FullSplice.Rng, mcols(splice.FullSplice.Rng)$matchRng, 
                         fill.gap=TRUE) == splice.FullSplice.Rng
  firstExonInIntron <- punion(mcols(start.Splice.Rng)$matchRng, start.Splice.Rng, 
                              fill.gap=TRUE) == mcols(start.Splice.Rng)$matchRng
  lastExonInIntron <- punion(mcols(start.Splice.Rng)$matchRng, end.Splice.Rng, 
                             fill.gap=TRUE) == mcols(start.Splice.Rng)$matchRng
  exonSkippingVector <- pmax(0, tapply(exonSkipping, mcols(splice.FullSplice.Rng)$IdMap, sum) -
                               tapply(firstExonInIntron, mcols(start.Splice.Rng)$IdMap, sum) - 
                               tapply(lastExonInIntron, mcols(start.Splice.Rng)$IdMap, sum))
  return (exonSkippingVector)
}

#' annotate exon splicing
#' @noRd
annotateExonSplice <- function(spliceRng, fullRng, startRng, endRng, strand){
  splice.FullSplice.Rng <- expandSpliceRanges(spliceRng, fullRng)
  start.Splice.Rng <- expandStartEndRanges(startRng, spliceRng)
  end.Splice.Rng <- rep(endRng, elementNROWS(spliceRng))
  fullSpliceTable <- tibble(splice.FullSplice.start = start(splice.FullSplice.Rng),
                            splice.FullSplice.end = end(splice.FullSplice.Rng),
                            match.FullSplice.start = start(mcols(splice.FullSplice.Rng)$matchRng),
                            match.FullSplice.end = end(mcols(splice.FullSplice.Rng)$matchRng))
  startEndSpliceTable <- tibble(startSplice.start = start(start.Splice.Rng),
                                startSplice.end = end(start.Splice.Rng),
                                match.startSplice.start = start(mcols(start.Splice.Rng)$matchRng),
                                match.startSplice.end = end(mcols(start.Splice.Rng)$matchRng),
                                endSplice.start = start(end.Splice.Rng),
                                endSplice.end = end(end.Splice.Rng))
  fullSpliceTable <- fullSpliceTable %>% mutate(exonStart = splice.FullSplice.end<match.FullSplice.end & 
                                                  splice.FullSplice.end>match.FullSplice.start &
                                                  splice.FullSplice.start<match.FullSplice.start,
                                                exonEnd = splice.FullSplice.start<match.FullSplice.end & 
                                                  splice.FullSplice.start>match.FullSplice.start &
                                                  splice.FullSplice.end>match.FullSplice.end)
  startEndSpliceTable <- startEndSpliceTable %>% mutate(startExonStartExtension = startSplice.start<match.startSplice.end &
                                                          startSplice.end>match.startSplice.end &
                                                          match.startSplice.start< startSplice.start,
                                                        startExonEndExtension = startSplice.end>match.startSplice.start &
                                                          startSplice.start<match.startSplice.start &
                                                          match.startSplice.end> startSplice.end,
                                                        endExonStartExtension = endSplice.start<match.startSplice.end &
                                                          endSplice.end>match.startSplice.end &
                                                          match.startSplice.start< endSplice.start,
                                                        endExonEndExtension = endSplice.end>match.startSplice.start &
                                                          endSplice.start<match.startSplice.start &
                                                          match.startSplice.end> endSplice.end)
  exon5Prime <- tapply(fullSpliceTable$exonStart, mcols(splice.FullSplice.Rng)$IdMap, sum)
  exon3Prime <- tapply(fullSpliceTable$exonEnd, mcols(splice.FullSplice.Rng)$IdMap, sum)
  startExonStartExtension <- tapply(startEndSpliceTable$startExonStartExtension, mcols(start.Splice.Rng)$IdMap, sum)
  startExonEndExtension <- tapply(startEndSpliceTable$startExonEndExtension, mcols(start.Splice.Rng)$IdMap, sum)
  endExonStartExtension <- tapply(startEndSpliceTable$endExonStartExtension, mcols(start.Splice.Rng)$IdMap, sum)
  endExonEndExtension <- tapply(startEndSpliceTable$endExonEndExtension, mcols(start.Splice.Rng)$IdMap, sum)
  exonSplicingTable <- tibble(exon5Prime,exon3Prime,strand)
  exStrandPos <- exonSplicingTable$strand!='-'
  exStrandNeg <- exonSplicingTable$strand=='-'
  exonSplicingTable$exon5Prime[exStrandPos] <- exonSplicingTable$exon5Prime[exStrandPos] - startExonStartExtension[exStrandPos]
  exonSplicingTable$exon5Prime[exStrandNeg] <- exon3Prime[exStrandNeg] - startExonEndExtension[exStrandNeg]
  exonSplicingTable$exon3Prime[exStrandPos] <- exonSplicingTable$exon3Prime[exStrandPos] - endExonEndExtension[exStrandPos]
  exonSplicingTable$exon3Prime[exStrandNeg] <- exon5Prime[exStrandNeg] - endExonStartExtension[exStrandNeg] 
  exonSplicingTable <- exonSplicingTable %>% select(exon5Prime,exon3Prime)
  return(exonSplicingTable)
}

#' Get intron ranges from exon ranges list
#' @noRd
unlistIntrons <- function(x, use.ids = TRUE, use.names=TRUE)
{
  # License note: This function is adopted from the GenomicAlignments package (Author: Hervé Pagès, Valerie Obenchain, Martin Morgan)
  # License	Artistic-2.0
  # https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments
  
  flat <- unlist(x, use.names=FALSE)
  gaps <- gaps(ranges(x))
  
  firstseg <- start(PartitioningByWidth(x))
  seqnms <- rep(seqnames(flat)[firstseg], elementNROWS(gaps))
  strand <- rep(strand(flat)[firstseg],elementNROWS(gaps))
 
  gr <- GRanges(seqnms, unlist(gaps, use.names=use.names), strand)
  if(use.ids & !is.null(mcols(x, use.names=FALSE)$id)) {
    mcols(gr)$id <- rep(mcols(x)$id,elementNROWS(gaps))
  }
  gr
}

#' Get intron ranges from exon ranges list
#' @noRd
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
#myGaps <- .GenomicAlignments:::.gaps

#' @noRd
.isNumericOrNAs <- S4Vectors:::isNumericOrNAs



#' @noRd
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
  DataFrame(uniqueLengthQuery,uniqueLengthSubject, compatible, queryElementsOutsideMaxDist, subjectElementsOutsideMaxDist)
}



#' The following function is implemented in R (GenomicAlignments), I just included the "within" option to make them significantly faster+memorey friendly for this purpose (original code modified from the GenomicAlignments package, Author: Hervé Pagès, Valerie Obenchain, Martin Morgan)
# License	Artistic-2.0
# https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments

#' @noRd
findSpliceOverlapsQuick <- function(query, subject, ignore.strand=FALSE) {
  olap <- findOverlaps(query, subject, ignore.strand=ignore.strand, type='within')
  olapEqual <- findOverlaps(query, subject, ignore.strand=ignore.strand, type='equal')
  if (length(olap) == 0L)
    return(GenomicAlignments:::.result(olap))


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

#' @param query query
#' @param subject subject
#' @param splice splice
#' @noRd
myCompatibleTranscription <- function(query, subject, splice)
{
  qrng <- ranges(query)
  srng <- ranges(subject)
  sprng <- ranges(splice)

  bnds <- elementNROWS(GenomicRanges::setdiff(qrng, srng)) == 0L
  splc <- elementNROWS(GenomicRanges::intersect(srng, sprng)) == 0L
  return(bnds & splc)
}

#' @param idx idx
#' @param x x
#' @noRd
myOneMatch <- function(x, idx)
{
  # License note: This function is adopted from the GenomicAlignments package (Author: Hervé Pagès, Valerie Obenchain, Martin Morgan)
  # https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments
  xcnt <- rowsum(as.integer(x), idx)[,1]
  oneMatch <- rep((xcnt == 1L), table(idx))
  unname(x & oneMatch)
}

#' @param motif motif
#' @noRd
spliceStrand <- function(motif){
  NATURAL_INTRON_MOTIFS_RC <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(GenomicAlignments::NATURAL_INTRON_MOTIFS)))

  motifStrand <- ifelse(motif %in% GenomicAlignments::NATURAL_INTRON_MOTIFS,'+','*')
  motifStrand[motif %in% NATURAL_INTRON_MOTIFS_RC] <- '-'
  return(motifStrand)
}


#' Function to reduce the start end end of the first and last elements in a granges list objects to a single basepair, helper to identify overlaps based on splicing only (allow for flexible TSS/TES)
#' @param grangesList grangesList
#' @noRd
cutStartEndFromGrangesList <- function(grangesList) {
  unlistedExons <- unlist(grangesList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(grangesList)), names=NULL)
  startExonsSet <- (which((unlistedExons$exon_rank==1 & as.character(strand(unlistedExons))!='-')|(unlistedExons$exon_endRank==1 & as.character(strand(unlistedExons))=='-')))
  endExonsSet <- (which((unlistedExons$exon_rank==1 & as.character(strand(unlistedExons))=='-')|(unlistedExons$exon_endRank==1 & as.character(strand(unlistedExons))!='-')))

  start(unlistedExons[startExonsSet]) <-  end(unlistedExons[startExonsSet])-1
  end(unlistedExons[endExonsSet]) <-  start(unlistedExons[endExonsSet])+1

  return(relist(unlistedExons, partitioning))
}

#' @param grangesList grangesList
#' @param by defaults to 5
#' @noRd
extendGrangesListElements <- function(grangesList, by=5) {
  unlistedExons <- unlist(grangesList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(grangesList)), names=NULL)
  start(unlistedExons) <- pmax(1,start(unlistedExons)-by)
  end(unlistedExons) <- pmin(seqlengths(unlistedExons)[as.character(seqnames(unlistedExons))],
                             end(unlistedExons)+by,
                             na.rm=TRUE)

  return(relist(unlistedExons, partitioning))
}


#' @param grangesList grangesList
#' @param minWidth defaults to 5
#' @param cutStartEnd defaults to FALSE
#' @noRd
dropGrangesListElementsByWidth <- function(grangesList, minWidth=5, cutStartEnd=FALSE) {
  unlistedExons <- unlist(grangesList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(sum(width(grangesList)>=minWidth)), names=NULL)
  exonWidth <- width(unlistedExons)
  if(cutStartEnd) {
    startExonsSet <- (which((unlistedExons$exon_rank==1 & as.character(strand(unlistedExons))!='-')|(unlistedExons$exon_endRank==1 & as.character(strand(unlistedExons))=='-')))
    endExonsSet <- (which((unlistedExons$exon_rank==1 & as.character(strand(unlistedExons))=='-')|(unlistedExons$exon_endRank==1 & as.character(strand(unlistedExons))!='-')))

    start(unlistedExons[startExonsSet]) <-  end(unlistedExons[startExonsSet])-1
    end(unlistedExons[endExonsSet]) <-  start(unlistedExons[endExonsSet])+1
  }
  unlistedExons <- unlistedExons[exonWidth>= minWidth]

  return(relist(unlistedExons, partitioning))
}

#' Function that selects the first exon from an IRangesList object
#' @param range IRangesList with elements required to be ordered by coordinates
#' @param stand strand
#' @noRd
selectStartExonFromRangesList <- function(range, strand){
  largestExons <- as.numeric(cumsum(elementNROWS(range)))
  startExonsSet <- c(1, largestExons[-(length(largestExons))]+1)
  startExonsSet[strand == "-"] <- largestExons[strand == "-"]
  return(unlist(range, use.names = FALSE)[startExonsSet])
}

#' Function that selects the last exon from an IRangesList object
#' @param range IRangesList with elements required to be ordered by coordinates
#' @param stand strand
#' @noRd
selectEndExonFromRangesList <- function(range, strand){
  endExonsSet <- as.numeric(cumsum(elementNROWS(range)))
  smallestExons <- c(1, endExonsSet[-(length(endExonsSet))]+1)
  endExonsSet[strand == "-"] <- smallestExons[strand == "-"]
  return(unlist(range, use.names = FALSE)[endExonsSet])
}
####the ones with relist()
# selectStartExonFromRangesList <- function(range, strand){
#   partitioning<-PartitioningByEnd(cumsum(pmin(elementNROWS(range),1)), names=NULL)
#   largestExons <- as.numeric(cumsum(elementNROWS(range)))
#   startExonsSet <- c(1, largestExons[1:length(largestExons)-1]+1)
#   startExonsSet[strand == "-"] <- largestExons[strand == "-"]
#   return(relist(unlist(range, use.names = FALSE)[startExonsSet],partitioning))
# }
# selectEndExonFromRangesList <- function(range, strand){
#   partitioning<-PartitioningByEnd(cumsum(pmin(elementNROWS(range),1)), names=NULL)
#   endExonsSet <- as.numeric(cumsum(elementNROWS(range)))
#   smallestExons <- c(1, endExonsSet[1:length(endExonsSet)-1]+1)
#   endExonsSet[strand == "-"] <- smallestExons[strand == "-"]
#   return(relist(unlist(range, use.names = FALSE)[endExonsSet],partitioning))
# }

#' Function that selects the first N exons from a grangeslist object (exon_rank is required)
#' @param grangesList grangesList
#' @param exonNumber defaults to 2
#' @noRd
selectStartExonsFromGrangesList <- function(grangesList, exonNumber=2) {
  unlisted_granges <- unlist(grangesList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(pmin(elementNROWS(grangesList), exonNumber)), names=NULL)
  startExonsSet <- which(unlisted_granges$exon_rank<=exonNumber)
  return(relist(unlisted_granges[startExonsSet], partitioning))
}

#' Function that selects the last N exons from a grangeslist object (exon_endRank is required)
#' @describeIn selectStartExonsFromGrangesList grangesList
#' @noRd
selectEndExonsFromGrangesList <- function(grangesList, exonNumber=2) {
  unlisted_granges <- unlist(grangesList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(pmin(elementNROWS(grangesList), exonNumber)), names=NULL)
  endExonsSet <- which(unlisted_granges$exon_endRank<=exonNumber)
  return(relist(unlisted_granges[endExonsSet], partitioning))
}


