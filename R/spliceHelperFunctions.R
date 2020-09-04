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
  queryStartRng <- ranges(selectStartExonsFromGrangesList(query, exonNumber = 1))
  subjectStartRng <- ranges(selectStartExonsFromGrangesList(subject, exonNumber = 1))
  queryEndRng<- ranges(selectEndExonsFromGrangesList(query, exonNumber = 1))
  subjectEndRng <- ranges(selectEndExonsFromGrangesList(subject, exonNumber = 1))
  subjectFullRng <- ranges(subject)
  queryFullRng <- ranges(query)
  querySpliceRng <- ranges(myGaps(query))
  querySpliceRng[elementNROWS(querySpliceRng)==0] <- IRanges(start=1,end=1) # add mock intron
  subjectSpliceRng <- ranges(myGaps(subject))
  subjectSpliceRng[elementNROWS(subjectSpliceRng)==0] <- IRanges(start=1,end=1)# add mock intron
  ## start end 
  annotatedTable <- tibble(queryId = names(query),
                           subjectId = names(subject),
                           strand = as.character(unlist(unique(strand(query)))))
  ###############################PRE-PROCESSING#################################
  #for intron retention/exon skipping
  querySplice.subjectFullQuerySplice.Rng <- expandMatchedRanges(querySpliceRng, subjectFullRng)
  subjectSplice.queryFullsubjectSplice.Rng <- expandMatchedRanges(subjectSpliceRng ,queryFullRng)
  #for exon skipping
  subjectStart.querySplice.Rng <- expandMatchedRanges(subjectStartRng, querySpliceRng)
  subjectEnd.querySplice.Rng <- rep(unlist(subjectEndRng), elementNROWS(querySpliceRng))
  queryStart.subjectSplice.Rng <- expandMatchedRanges(queryStartRng, subjectSpliceRng)
  queryEnd.subjectSplice.Rng <- rep(unlist(queryEndRng), elementNROWS(subjectSpliceRng))
  ###############################################################################
  # calculate alternative First/last exons and annotate internal start and end first exons
  alterEndTable <- annotateAlterEndExons(queryStartRng, queryEndRng, queryFullRng,
                                         subjectStartRng, subjectEndRng, subjectFullRng, annotatedTable$strand)
  #intron retention subject
  annotatedTable$intronRetention.subject <- annotateIntronRetent(querySplice.subjectFullQuerySplice.Rng)
  #intron retention query
  annotatedTable$intronRetention.query <- annotateIntronRetent(subjectSplice.queryFullsubjectSplice.Rng)
  #exon skipping query
  annotatedTable$exonSkipping.query <- annotateExonSkip(querySplice.subjectFullQuerySplice.Rng,
                                                        subjectStart.querySplice.Rng,
                                                        subjectEnd.querySplice.Rng)
  #exon skipping subject
  annotatedTable$exonSkipping.subject <- annotateExonSkip(subjectSplice.queryFullsubjectSplice.Rng,
                                                          queryStart.subjectSplice.Rng,
                                                          queryEnd.subjectSplice.Rng) 
  # exon 3' splice site
  exonSpliceTable <- annotateExonSplice(querySplice.subjectFullQuerySplice.Rng,
                                        subjectStart.querySplice.Rng,
                                        subjectEnd.querySplice.Rng,
                                        annotatedTable$strand)
  annotatedTable <- cbind(annotatedTable, alterEndTable, exonSpliceTable)
  return(annotatedTable)
}

#' ranges pre-processing                                  
expandMatchedRanges <- function(Rng,comparedRng){ 
  FinRng <- rep(unlist(Rng, use.names=F),rep(elementNROWS(comparedRng),times=elementNROWS(Rng)))
  mcols(FinRng)$IdMap <- rep(1:length(Rng),elementNROWS(Rng)*elementNROWS(comparedRng))
  mcols(FinRng)$matchRng <- unlist(rep(comparedRng, times= elementNROWS(Rng)), use.names=F)
  return (FinRng)
}

#' annotate alternative First/last exons and TSS/TES                                  
#' @noRd
annotateAlterEndExons <- function(queryStartRng, queryEndRng, queryFullRng,
                                  subjectStartRng, subjectEndRng, subjectFullRng, strand){
  tempTable <-  tibble(start.first.query = start(unlist(queryStartRng)),
                       end.first.query = end(unlist(queryStartRng)),
                       start.last.query = start(unlist(queryEndRng)),
                       end.last.query = end(unlist(queryEndRng)),
                       start.first.subject = start(unlist(subjectStartRng)),
                       end.first.subject = end(unlist(subjectStartRng)),
                       start.last.subject = start(unlist(subjectEndRng)),
                       end.last.subject = end(unlist(subjectEndRng)),
                       strand = strand)  
  tempTable <- tempTable %>% mutate(alternativeFirstExon=!ifelse(strand!='-',
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
  tempTable <- tempTable %>% select(alternativeFirstExon, alternativeLastExon, alternativeTSS, alternativeTES)
  #internal start/end query 
  queryTable <- annotateInternalStartEnd(queryStartRng, queryEndRng, subjectFullRng, tempTable)
  colnames(queryTable) <- paste(colnames(queryTable),"query",sep=".")
  #internal start/end subject
  subjectTable <- annotateInternalStartEnd(subjectStartRng, subjectEndRng, queryFullRng, tempTable)
  colnames(subjectTable) <- paste(colnames(subjectTable),"subject",sep=".")
  tempTable <- cbind(tempTable, queryTable, subjectTable)
  return (tempTable)
}

#' annotate internal start and end first exons                                    
#' @noRd
annotateInternalStartEnd <- function(StartRng, EndRng, FullRng, alterEndTable){
  Start.Full.Rng <- expandMatchedRanges(StartRng, FullRng) 
  End.Full.Rng <- expandMatchedRanges(EndRng, FullRng)
  #internal start
  ExonIntersectStart <- pintersect(Start.Full.Rng,mcols(Start.Full.Rng)$matchRng, resolve.empty='start.x')
  #internal end
  ExonIntersectEnd=pintersect(End.Full.Rng,mcols(End.Full.Rng)$matchRng, resolve.empty='start.x')
  tempTable <- tibble(internalFirstExon = tapply(width(ExonIntersectStart), mcols(Start.Full.Rng)$IdMap, sum)!= 0 & alterEndTable$alternativeFirstExon,
                      internalLastExon = tapply(width(ExonIntersectEnd), mcols(End.Full.Rng)$IdMap, sum)!=0 & alterEndTable$alternativeLastExon)
  return(tempTable)
}

#' annotate intron retention 
#' @noRd
annotateIntronRetent <- function(Splice.FullSplice.Rng, annotatedTable){
  IntronRetention <-punion(Splice.FullSplice.Rng, mcols(Splice.FullSplice.Rng)$matchRng, fill.gap=TRUE) == mcols(Splice.FullSplice.Rng)$matchRng
  outputVector <- tapply(IntronRetention, mcols(Splice.FullSplice.Rng)$IdMap, sum)
  return (outputVector)
}

#' annotate exon skiping 
#' @noRd
annotateExonSkip <- function(Splice.FullSplice.Rng, Start.Splice.Rng, End.Splice.Rng){
  ExonSkipping <- punion(Splice.FullSplice.Rng, mcols(Splice.FullSplice.Rng)$matchRng, 
                         fill.gap=TRUE) == Splice.FullSplice.Rng
  FirstExonInIntron <- punion(mcols(Start.Splice.Rng)$matchRng, Start.Splice.Rng, 
                              fill.gap=TRUE) == mcols(Start.Splice.Rng)$matchRng
  LastExonInIntron <- punion(mcols(Start.Splice.Rng)$matchRng, End.Splice.Rng, 
                             fill.gap=TRUE) == mcols(Start.Splice.Rng)$matchRng
  outputVector <- pmax(0, tapply(ExonSkipping, mcols(Splice.FullSplice.Rng)$IdMap, sum) -
                          tapply(FirstExonInIntron, mcols(Start.Splice.Rng)$IdMap, sum) - 
                          tapply(LastExonInIntron, mcols(Start.Splice.Rng)$IdMap, sum))
  return (outputVector)
}

#' annotate exon splicing
#' @noRd
annotateExonSplice <- function(Splice.FullSplice.Rng, Start.Splice.Rng,End.Splice.Rng, strand){
  ExonStart <- end(Splice.FullSplice.Rng)<end(mcols(Splice.FullSplice.Rng)$matchRng) & 
    end(Splice.FullSplice.Rng)>start(mcols(Splice.FullSplice.Rng)$matchRng) &
    start(Splice.FullSplice.Rng)<start(mcols(Splice.FullSplice.Rng)$matchRng)
  ExonEnd <- start(Splice.FullSplice.Rng)<end(mcols(Splice.FullSplice.Rng)$matchRng) & 
    start(Splice.FullSplice.Rng)>start(mcols(Splice.FullSplice.Rng)$matchRng) &
    end(Splice.FullSplice.Rng)>end(mcols(Splice.FullSplice.Rng)$matchRng)
  
  exon5Prime <- tapply(ExonStart, mcols(Splice.FullSplice.Rng)$IdMap, sum)
  exon3Prime <- tapply(ExonEnd, mcols(Splice.FullSplice.Rng)$IdMap, sum)
  
  StartExonStartExtension <- start(Start.Splice.Rng)<end(mcols(Start.Splice.Rng)$matchRng) &
    end(Start.Splice.Rng)>end(mcols(Start.Splice.Rng)$matchRng) &
    start(mcols(Start.Splice.Rng)$matchRng)< start(Start.Splice.Rng)
  StartExonStartExtension <- tapply(StartExonStartExtension, mcols(Start.Splice.Rng)$IdMap, sum)
  
  StartExonEndExtension <- end(Start.Splice.Rng)>start(mcols(Start.Splice.Rng)$matchRng) &
    start(Start.Splice.Rng)<start(mcols(Start.Splice.Rng)$matchRng) &
    end(mcols(Start.Splice.Rng)$matchRng)> end(Start.Splice.Rng)
  StartExonEndExtension <- tapply(StartExonEndExtension, mcols(Start.Splice.Rng)$IdMap, sum)
  
  EndExonStartExtension <- start(End.Splice.Rng)<end(mcols(Start.Splice.Rng)$matchRng) &
    end(End.Splice.Rng)>end(mcols(Start.Splice.Rng)$matchRng) &
    start(mcols(Start.Splice.Rng)$matchRng)< start(End.Splice.Rng)
  EndExonStartExtension <- tapply(EndExonStartExtension, mcols(Start.Splice.Rng)$IdMap, sum)
  
  EndExonEndExtension <- end(End.Splice.Rng)>start(mcols(Start.Splice.Rng)$matchRng) &
    start(End.Splice.Rng)<start(mcols(Start.Splice.Rng)$matchRng) &
    end(mcols(Start.Splice.Rng)$matchRng)> end(End.Splice.Rng)
  EndExonEndExtension <- tapply(EndExonEndExtension, mcols(Start.Splice.Rng)$IdMap, sum)
  
  tempTable <- tibble(exon5Prime,exon3Prime,strand)
  tempTable$exon5Prime[tempTable$strand!='-'] <- tempTable$exon5Prime[tempTable$strand!='-'] - StartExonStartExtension[tempTable$strand!='-']
  tempTable$exon5Prime[tempTable$strand=='-'] <- exon3Prime[tempTable$strand=='-'] - StartExonEndExtension[tempTable$strand=='-']

  tempTable$exon3Prime[tempTable$strand!='-'] <- tempTable$exon3Prime[tempTable$strand!='-'] - EndExonEndExtension[tempTable$strand!='-']
  tempTable$exon3Prime[tempTable$strand=='-'] <- exon5Prime[tempTable$strand=='-'] - EndExonStartExtension[tempTable$strand=='-'] 
  tempTable <- tempTable %>% select(exon5Prime,exon3Prime)
  return(tempTable)
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


