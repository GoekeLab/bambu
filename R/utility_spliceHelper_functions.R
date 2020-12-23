# Note: several of the functions in this file are adopted from 
# the GenomicAlignments package (Author: Hervé Pagès, Valerie Obenchain,
# Martin Morgan)
# License Artistic-2.0
# https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments

### examples for test purposes
## Expected annotations of transcripts used in test query
# 'ENST00000344579', # exon skipping, alternative TSS (-48), +, ENSG00000158109
# 'ENST00000270792', # intron retention subject 1(last exon),alt.TSS,alt.TES, +,
# 'ENST00000410032', # alternative first exon, exon skipping query: 2, 
#                    # exon skipping subject: 0, alternative TSS (2bp only), 
#                    # internalFirstExon.subject +
# 'ENST00000468178', # alternative last exon +
# 'ENST00000485956', # alternative first exon, alternative last exon,
# #exon skipping subject = 1, internal first exon query, +
# 'ENST00000530807', # exon skipping query 1, alternative TSS (-17),  -
# 'ENST00000409894', # alternative 3' exon splice site, exon skipping query 2,
# #alternative TSS, alterantive TES, +, ENSG00000125630
# 'ENST00000524447',  # alternative TSS, alternative last exon (internal), 
# #alternative exon 3' end,-, ENSG00000165916
# 'ENST00000591696' # alternative TSS, alternative 3' exon (2), 
# #alternative 5' exon (1) alternative TES, ,+,ENSG00000141349
############################################################
# query <- readRDS(system.file("extdata", 
#     "annotateSpliceOverlapByDist_testQuery.rds",
#     package = "bambu"))
# subject <- readRDS(system.file("extdata", 
#        "annotateSpliceOverlapByDist_testSubject.rds",
#        package = "bambu"))
############################################################

#' annotate splice overlap by distance
#' @description This function takes in a GRangesList (query)
#' and a target GRangesList (subject). The function creates
#' an annotation table in tibble by comparing ranges entries
#' from transcripts between the query and subject GRangesLists.
#' @usage compareTranscripts(query, subject)
#' @params query a GRangesList
#' @params subject a GRangesList
#' @return a tibble with the following annotations:
#' \itemize{
#'    \item alternativeFirstExon
#'    \item alternativeTSS
#'    \item internalFirstExon.query
#'    \item internalFirstExon.subject
#'    \item alternativeLastExon
#'    \item alternativeTES
#'    \item internalLastExon.query
#'    \item internalLastExon.subject
#'    \item intronRetention.subject
#'    \item intronRetention.query
#'    \item exonSkipping.query
#'    \item exonSkipping.subject
#'    \item exon5prime (splicing)
#'    \item exon3prime (splicing)
#' }
#' @examples
#' query <- readRDS(system.file("extdata", 
#'     "annotateSpliceOverlapByDist_testQuery.rds",
#'     package = "bambu"))
#' subject <- readRDS(system.file("extdata", 
#'     "annotateSpliceOverlapByDist_testSubject.rds",
#'     package = "bambu"))
#' annotationTable <- compareTranscripts(query, subject)
#' @noRd
compareTranscripts <-function(query, subject) {
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
  annotatedTable <- tibble(queryId = names(query), subjectId = names(subject), strand = strand)
  # calculate alternative First/last exons and annotate internal start and end first exons
  annotatedTable$alternativeFirstExon <- alternativeStartEndExon(queryStartRng, 
                                                                 subjectStartRng)
  annotatedTable$alternativeTSS <- calculateTerminalDistance(queryStartRng, 
                                                        subjectStartRng, 
                                                        annotatedTable$alternativeFirstExon,
                                                        strand, "start")
  annotatedTable$internalFirstExon.query <- annotateInternalStartEnd(queryStartRng,
                                                                     subjectFullRng, annotatedTable$alternativeFirstExon)
  annotatedTable$internalFirstExon.subject <- annotateInternalStartEnd(subjectStartRng,
                                                                       queryFullRng, annotatedTable$alternativeFirstExon)
  # TES/last exon
  annotatedTable$alternativeLastExon <- alternativeStartEndExon(queryEndRng, 
                                                                subjectEndRng)
  annotatedTable$alternativeTES <- calculateTerminalDistance(queryEndRng, 
                                                        subjectEndRng, annotatedTable$alternativeLastExon, strand, "end")
  annotatedTable$internalLastExon.query <- annotateInternalStartEnd(queryEndRng,
                                                                    subjectFullRng, annotatedTable$alternativeLastExon)
  annotatedTable$internalLastExon.subject <- annotateInternalStartEnd(subjectEndRng,
                                                                      queryFullRng, annotatedTable$alternativeLastExon)
  #intron retention subject and query
  annotatedTable$intronRetention.subject <- annotateIntronRetent(querySpliceRng,
                                                                 subjectFullRng)
  annotatedTable$intronRetention.query <- annotateIntronRetent(subjectSpliceRng,
                                                               queryFullRng)
  #exon skipping query and subject
  annotatedTable$exonSkipping.query <- annotateExonSkip(querySpliceRng, subjectFullRng,
                                                        subjectStartRng, subjectEndRng)
  annotatedTable$exonSkipping.subject <- annotateExonSkip(subjectSpliceRng ,queryFullRng,
                                                          queryStartRng, queryEndRng) 
  #exon 5' and 3' splice site
  exonSpliceTable <- annotateExonSplice(querySpliceRng, subjectFullRng,
                                        subjectStartRng, subjectEndRng,
                                        annotatedTable$strand)
  annotatedTable <- cbind(annotatedTable, exonSpliceTable)
  return(annotatedTable)
}

#' extract strand from GRangesList
#' @description This function takes a GRangesList and
#' returns a vector with the strand for each list entry.
#' This function assumes that all elements for each list
#' entry have the same strand.
#' @usage getStrangeFromGrList(grl)
#' @params grl a GRangesList
#' @return an Rle object with strand information
#' @examples
#' query <- readRDS(system.file("extdata", 
#'     "annotateSpliceOverlapByDist_testQuery.rds",
#'     package = "bambu"))
#' strand <- as.character(getStrandFromGrList(query))
#' @noRd
getStrandFromGrList <- function(grl) { 
  return(unlist(strand(grl), use.names = F)[cumsum(elementNROWS(grl))]) 
}

#' start/end ranges pre-processing
#' @description this function takes in an IRanges object and a target
#' IRangesList object with the same length, where each list entry
#' i in IRangesList (target[[i]]) corresponds to the matching range
#' i in IRanges (ranges[i]). The function then creates a new IRanges
#' object with a length corresponding to length(unlist(target)), 
#' where the ranges[i] elements are repeated to match each individual
#' element in target[[i]]. The unlist(target) ranges are stored in
#' mcols()$matchRng while the corresponding index i for each element
#' of target[[i]] is stored in mcols()$IdMap. This function is
#' used to enable the comparison of ranges with all elements in a
#' rangesList, for example to compute the overlap of first and last
#' exons with matching transcripts (see compareTranscripts()).
#' @param ranges an IRanges object
#' @param target an IRangesList object
#' @return a ranges object with mcols objects
#' \itemize{
#'    \item matchRng repeated matched ranges
#'    \item IdMap index of the repeated matched ranges
#' }
#' @noRd
expandRanges <- function(ranges,target){ 
  processedRng <- rep(ranges,elementNROWS(target))
  mcols(processedRng)$IdMap <- rep(1:length(ranges),elementNROWS(target))
  mcols(processedRng)$matchRng <- unlist(target, use.names=F)
  return (processedRng)
}

#' splice ranges pre-processing
#' @description this function takes in an IRangesList object and a 
#' target IRangesList object with the same length, where each list
#' entry in target IRangesList[[i]] (target[[i]]) corresponds to the
#' matching list entry in IRangesList[[i]] (rglist[i]). The function
#' then creates an IRanges object with a length corresponding to the
#' length(unlist(target))*length(unlist(rglist)), where each individual
#' element in rglist[[i]] is repeated to match each individual element
#' in target[[i]]. The repeated ranges are stored in mcols()$matchRng 
#' while the corresponding index i for each element of 
#' target[[i]]*rglist[[i]] is stored in mcols()$IdMap. This function is
#' used to enable the comparison of ranges with all elements in a
#' rangesList, for example to compute the overlap of splice sites
#' with matching transcripts (see compareTranscripts()).
#' @param rglist an IRangesList object
#' @param target an IRangesList object
#' @return a ranges object with mcols objects
#' \itemize{
#'    \item matchRng repeated matched ranges
#'    \item IdMap index of the repeated matched ranges
#' }
#' @noRd
expandRangesList <- function(rglist,target){ 
  processedRng <- rep(unlist(rglist, use.names=F),
                      rep(elementNROWS(target),times=elementNROWS(rglist)))
  mcols(processedRng)$IdMap <- rep(1:length(rglist),
                                   elementNROWS(rglist)*elementNROWS(target))
  mcols(processedRng)$matchRng <- unlist(rep(target, 
                                             times= elementNROWS(rglist)), use.names=F)
  return (processedRng)
}

#' alternative start/end exon
#' @description This function checks whether an 
#' alternative start/end exon is used by overlapping
#' the exon ranges of the first (or last) exons of matching transcripts.
#' @noRd
alternativeStartEndExon <- function(queryRng, subjectRng){
  return(!poverlaps(queryRng, subjectRng))
}

#' alternative TSS/TES distance
#' @description This function calculates the distance of an
#' alternative TSS/TES by comparing the start/end coordinates
#' of the start/end exon ranges of matching transcripts. If an 
#' alternative first/last exon is used the distance is set to 0.
#' @noRd
calculateTerminalDistance <- function(queryTerminalExonRng,
                                      subjectTerminalExonRng, alternativeTerminalExon,
                                      strand, direction = "start"){
    direction_names <- c("start","end")
    alternativeTerminal <- (-1)^(direction == "end") *
           (get(direction)(subjectTerminalExonRng) - get(direction)(queryTerminalExonRng))
    alternativeTerminal[strand == "-"] <- (-1)^(direction == "end") * 
      (get(setdiff(direction_names, direction))(queryTerminalExonRng[strand=='-'])-
       get(setdiff(direction_names, direction))(subjectTerminalExonRng[strand=='-']))
    alternativeTerminal <- alternativeTerminal * !alternativeTerminalExon
    return(alternativeTerminal)
}

#' annotate internal start and end first exons
#' @description This function checks whether
#' there is an internal start/end by overlapping
#' the exon ranges of matching transcripts.
#' @noRd
annotateInternalStartEnd <- function(exonRng, fullRng, 
                                     alternativeFirstLastExon){
  exon.Full.Rng <- expandRanges(exonRng, fullRng) 
  #internal start/end
  exonIntersect <- pintersect(exon.Full.Rng,
                              mcols(exon.Full.Rng)$matchRng, resolve.empty='start.x')
  internalStartEndVector <- tapply(width(exonIntersect),
                                   mcols(exon.Full.Rng)$IdMap, sum)!= 0 & alternativeFirstLastExon
  return(internalStartEndVector)
}

#' annotate intron retention
#' @description This function checks whether
#' there is intron retention by overlapping
#' the intron ranges of matching transcripts.
#' @noRd
annotateIntronRetent <- function(spliceRng, fullRng){
  splice.FullSplice.Rng <- expandRangesList(spliceRng, fullRng)
  intronRetention <-punion(splice.FullSplice.Rng, 
                           mcols(splice.FullSplice.Rng)$matchRng, fill.gap=TRUE) ==
    mcols(splice.FullSplice.Rng)$matchRng
  intronRetentionVector <- tapply(intronRetention,
                                  mcols(splice.FullSplice.Rng)$IdMap, sum)
  return (intronRetentionVector)
}

#' annotate exon skiping
#' @description This function checks whether
#' there is exon skipping by overlapping 
#' the intron ranges of matching transcripts.
#' @noRd
annotateExonSkip <- function(spliceRng, fullRng, startRng, endRng){
  splice.FullSplice.Rng <- expandRangesList(spliceRng, fullRng)
  start.Splice.Rng <- expandRanges(startRng, spliceRng)
  end.Splice.Rng <- rep(endRng, elementNROWS(spliceRng))
  exonSkipping <- punion(splice.FullSplice.Rng, 
                         mcols(splice.FullSplice.Rng)$matchRng, 
                         fill.gap=TRUE) == splice.FullSplice.Rng
  firstExonInIntron <- punion(mcols(start.Splice.Rng)$matchRng, start.Splice.Rng, 
                              fill.gap=TRUE) == mcols(start.Splice.Rng)$matchRng
  lastExonInIntron <- punion(mcols(start.Splice.Rng)$matchRng, end.Splice.Rng, 
                             fill.gap=TRUE) == mcols(start.Splice.Rng)$matchRng
  exonSkippingVector <- pmax(0, tapply(exonSkipping,
                                       mcols(splice.FullSplice.Rng)$IdMap, sum) -
                               tapply(firstExonInIntron, 
                                      mcols(start.Splice.Rng)$IdMap, sum) - 
                               tapply(lastExonInIntron, 
                                      mcols(start.Splice.Rng)$IdMap, sum))
  return (exonSkippingVector)
}

#' find exon start extension
#' @description This function checks whether
#' there is an extension at the start of an exon
#' by comparing the coordinates of the splice sites.
#' @noRd
findExonStartExtension <- function(splice.Rng, match.startSplice.start,
                                   match.startSplice.end, spliceIdMap){
  splice.start <- start(splice.Rng)
  splice.end <- end(splice.Rng)
  exonStartExtension <- splice.start<match.startSplice.end &
    splice.end>match.startSplice.end &
    match.startSplice.start< splice.start
  exonStartExtension <- tapply(exonStartExtension, spliceIdMap, sum)
  return (exonStartExtension)
}

#' find exon end extension
#' @description This function checks whether
#' there is an extension at the end of an exon
#' by comparing the coordinates of the splice sites.
#' @noRd
findExonEndExtension <- function(splice.Rng, match.startSplice.start,
                                 match.startSplice.end, spliceIdMap){
  splice.start <- start(splice.Rng)
  splice.end <- end(splice.Rng)
  exonEndExtension <- splice.end>match.startSplice.start &
    splice.start<match.startSplice.start &
    match.startSplice.end> splice.end
  exonEndExtension <- tapply(exonEndExtension, spliceIdMap, sum)
  return (exonEndExtension)
}

#' annotate exon splicing
#' @description This function checks whether 
#' there is alternative splicing in the 5'/3'
#' end of an exon.
#' @noRd
annotateExonSplice <- function(spliceRng, fullRng, startRng, endRng, strand){
  splice.FullSplice.Rng <- expandRangesList(spliceRng, fullRng)
  start.Splice.Rng <- expandRanges(startRng, spliceRng)
  end.Splice.Rng <- rep(endRng, elementNROWS(spliceRng))
  startMatch <- poverlaps(start(splice.FullSplice.Rng),
                          mcols(splice.FullSplice.Rng)$matchRng) 
  endMatch <- poverlaps(end(splice.FullSplice.Rng),
                        mcols(splice.FullSplice.Rng)$matchRng)
  exonStart <- endMatch & !startMatch
  exonEnd <-startMatch & !endMatch
  exon5Prime <- tapply(exonStart, mcols(splice.FullSplice.Rng)$IdMap, sum)
  exon3Prime <- tapply(exonEnd, mcols(splice.FullSplice.Rng)$IdMap, sum)
  match.startSplice.start <- start(mcols(start.Splice.Rng)$matchRng)
  match.startSplice.end <- end(mcols(start.Splice.Rng)$matchRng)
  spliceIdMap <- mcols(start.Splice.Rng)$IdMap
  startExonStartExtension <- findExonStartExtension(start.Splice.Rng,
                                                    match.startSplice.start,
                                                    match.startSplice.end,
                                                    spliceIdMap)
  startExonEndExtension <- findExonEndExtension(start.Splice.Rng,
                                                match.startSplice.start,
                                                match.startSplice.end,
                                                spliceIdMap)
  endExonStartExtension <- findExonStartExtension(end.Splice.Rng,
                                                  match.startSplice.start,
                                                  match.startSplice.end,
                                                  spliceIdMap)
  endExonEndExtension <- findExonEndExtension(end.Splice.Rng,
                                              match.startSplice.start,
                                              match.startSplice.end,
                                              spliceIdMap)
  exonSplicingTable <- tibble(exon5Prime,exon3Prime,strand)
  exStrandNeg <- exonSplicingTable$strand=='-'
  exonSplicingTable$exon5Prime[!exStrandNeg] <- 
    exonSplicingTable$exon5Prime[!exStrandNeg] - 
    startExonStartExtension[!exStrandNeg]
  exonSplicingTable$exon5Prime[exStrandNeg] <- exon3Prime[exStrandNeg] - 
    startExonEndExtension[exStrandNeg]
  exonSplicingTable$exon3Prime[!exStrandNeg] <- 
    exonSplicingTable$exon3Prime[!exStrandNeg] - 
    endExonEndExtension[!exStrandNeg]
  exonSplicingTable$exon3Prime[exStrandNeg] <- exon5Prime[exStrandNeg] - 
    endExonStartExtension[exStrandNeg] 
  exonSplicingTable <- exonSplicingTable %>% select(exon5Prime,exon3Prime)
  return(exonSplicingTable)
}

#' calculate distance between first and last exon matches
#' @param queryExon a query start or end exon ranges
#' @param subjectExon a subject start or end exon ranges
#' @param subjectFull a full subject ranges object
#' @param subjectList a full subject list
#' @noRd
calculateFirstLastExonsDist <- function(queryExon, subjectExon,
                                        subjectFull, subjectList) {
  # distances corresponds to length of unique sequences in all matching exons
  # distance = width(element in query) + width(all matching elements 
  # in subject) - 2x width(intersection(query, subject)
  ## subjectList <- unlist(subjectFull)
  ExonMatch <- poverlaps(unlist(queryExon), unlist(subjectExon))
  queryExonList <- rep(unlist(queryExon), elementNROWS(subjectFull))
  myId <- rep(seq_along(queryExon), elementNROWS(subjectFull))
  byExonIntersect <- pintersect(queryExonList, subjectList,
                                resolve.empty = "start.x")
  ExonDist <- as.integer(tapply(width(subjectList) *
                                  (width(byExonIntersect) > 0) - 2 * width(byExonIntersect),
                                myId, sum) + width(unlist(queryExon)))
  uniqueExonLengthQuery <-
    sum(width(GenomicRanges::setdiff(queryExon, subjectFull)))
  uniqueExonLengthSubject <- ExonDist - uniqueExonLengthQuery
  ExonMatchList <- list(
    "match" = ExonMatch,
    "uniqueExonLengthQuery" = uniqueExonLengthQuery,
    "uniqueExonLengthSubject" = uniqueExonLengthSubject)
  return(ExonMatchList)
}


#' This function calcualtes compatible splice overlaps allowing for a 
#' distance threshold, and returns distance in bp between query and subject.
#' Can be used to assign more transcripts to annotations and reads to
#' transcripts.
#' @noRd
findSpliceOverlapsByDist <- function(query, subject, ignore.strand = FALSE,
                                     maxDist = 5, type = "within", firstLastSeparate = TRUE,
                                     dropRangesByMinLength = FALSE, cutStartEnd = TRUE) {
  if (firstLastSeparate) {
    queryStart <- selectStartExonsFromGrangesList(query, exonNumber = 1)
    queryEnd <- selectEndExonsFromGrangesList(query, exonNumber = 1)
    subjectStart <- selectStartExonsFromGrangesList(subject, exonNumber = 1)
    subjectEnd <- selectEndExonsFromGrangesList(subject, exonNumber = 1)
    subjectFull <- subject
  }
  if (dropRangesByMinLength) {
    queryForOverlap <- dropGrangesListElementsByWidth(query,
                                                      minWidth = maxDist, cutStartEnd = cutStartEnd)
  } else if (cutStartEnd) {
    queryForOverlap <- cutStartEndFromGrangesList(query)
  } else {
    queryForOverlap <- query
  }
  query <- cutStartEndFromGrangesList(query)
  subjectExtend <- extendGrangesListElements(subject, by = maxDist)
  olap <- findOverlaps(queryForOverlap, subjectExtend,
                       ignore.strand = ignore.strand, type = type)
  olapEqual <- findOverlaps(query, cutStartEndFromGrangesList(subject),
                            ignore.strand = ignore.strand, type = "equal")
  query <- query[queryHits(olap)]
  subject <- subject[subjectHits(olap)]
  splice <- myGaps(query)
  compatible <- rangesDist(query, subject, splice, maxDist)
  equal <- (!is.na(S4Vectors::match(olap, olapEqual)))
  unique <- myOneMatch(compatible$compatible, queryHits(olap))
  strandSpecific <- all(strand(query) != "*")
  strandedMatch <- ((all(strand(query) == "-") & 
                       all(strand(subject) == "-")) | 
                      (all(strand(query) == "+") & 
                         all(strand(subject) == "+")))
  mcols(olap) <- DataFrame(compatible, equal, unique,
                           strandSpecific, strandedMatch)
  
  ## NOTE: Check if there is an error with the start sequence ##
  if (firstLastSeparate)
    olap <- checkStartSequence(olap, firstLastSeparate, queryStart,
                               subjectStart, queryEnd,subjectEnd, subjectFull, subjectList)
  return(olap)
}


#' check whether error with start sequence
#' @noRd
checkStartSequence <- function(olap, firstLastSeparate, queryStart,
                               subjectStart, queryEnd,subjectEnd, subjectFull, subjectList){
  if (length(olap)) {
    queryStart <- ranges(queryStart[queryHits(olap)])
    subjectStart <- ranges(subjectStart[subjectHits(olap)])
    queryEnd <- ranges(queryEnd[queryHits(olap)])
    subjectEnd <- ranges(subjectEnd[subjectHits(olap)])
    subjectFull <- ranges(subjectFull[subjectHits(olap)])
    subjectList <- unlist(subjectFull)
    startList <- calculateFirstLastExonsDist(queryStart, subjectStart,
                                             subjectFull, subjectList)
    endList <- calculateFirstLastExonsDist(queryEnd, subjectEnd,
                                           subjectFull, subjectList)
  } else {
    startList <- NULL
    endList <- NULL
  }
  mcols(olap) <- DataFrame(mcols(olap),
                           startMatch = startList$match,
                           uniqueStartLengthQuery = startList$uniqueExonLengthQuery,
                           uniqueStartLengthSubject = startList$uniqueExonLengthSubject,
                           endMatch = endList$match,
                           uniqueEndLengthQuery = endList$uniqueExonLengthQuery,
                           uniqueEndLengthSubject = endList$uniqueExonLengthSubject)
  return(olap)
}

#' Get intron ranges from exon ranges list
#' @noRd
unlistIntrons <- function(x, use.ids = TRUE, use.names = TRUE) {
  # License note: This function is adopted from the GenomicAlignments 
  # package (Author: Hervé Pagès, Valerie Obenchain, Martin Morgan)
  # License Artistic-2.0
  # https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments
  
  flat <- unlist(x, use.names = FALSE)
  gaps <- gaps(ranges(x))
  
  firstseg <- start(PartitioningByWidth(x))
  seqnms <- rep(seqnames(flat)[firstseg], elementNROWS(gaps))
  strand <- rep(strand(flat)[firstseg], elementNROWS(gaps))
  
  gr <- GenomicRanges::GRanges(seqnms, unlist(gaps,
                                              use.names = use.names), strand)
  if (use.ids & !is.null(mcols(x, use.names = FALSE)$id)) 
    mcols(gr)$id <- rep(mcols(x)$id, elementNROWS(gaps))
  return(gr)
}

#' Get intron ranges from exon ranges list
#' @noRd
myGaps <- function(x, start = NA, end = NA) {
  # License note: This function is adopted from the GenomicAlignments package
  # (Author: Hervé Pagès, Valerie Obenchain, Martin Morgan)
  # License Artistic-2.0
  # https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments
  if (!.isNumericOrNAs(start)) stop("'start' must be an integer vector or NA")
  if (!is.integer(start)) start <- as.integer(start)
  if (!.isNumericOrNAs(end)) stop("'end' must be an integer vector or NA")
  if (!is.integer(end)) end <- as.integer(end)
  
  ## seqname and strand consistent in list elements
  if (all(elementNROWS(runValue(seqnames(x))) == 1L) &&
      all(elementNROWS(runValue(strand(x))) == 1L)) {
    flat <- unlist(x, use.names = FALSE)
    gaps <- gaps(ranges(x), start, end)
    ### FIXME: this makes this function more of an 'introns' than a .gaps.
    ### FIXME: this breaks when the GRangesList is not ordered by position
    if (!is.null(mcols(x, use.names = FALSE)$query.break)) {
      insert_gaps <-
        methods::as(ranges(.insertGaps(x)), "CompressedIRangesList")
      gaps <- setdiff(gaps, insert_gaps)
    }
    
    idx <- elementNROWS(gaps) != 0
    ## FIXME : can't handle lists with empty elements
    ##         'start' and 'end' not quite right here
    firstseg <- start(PartitioningByWidth(x))
    seqnms <- rep(seqnames(flat)[firstseg], elementNROWS(gaps))
    strand <- rep(strand(flat)[firstseg], elementNROWS(gaps))
    gr <- relist(GenomicRanges::GRanges(seqnms, unlist(gaps,
                                                       use.names = FALSE), strand), gaps)
    gr
  } else {
    ### FIXME: does not handle query.break column yet
    setdiff(range(x), x)
  }
}
# myGaps <- .GenomicAlignments:::.gaps

#' @noRd
.isNumericOrNAs <- S4Vectors:::isNumericOrNAs



#' @noRd
rangesDist <- function(query, subject, splice, maxDist) {
  qrng <- ranges(query)
  srng <- ranges(subject)
  sprng <- ranges(splice)
  
  setDiffQ <- width(GenomicRanges::setdiff(qrng, srng))
  interesectS <- width(GenomicRanges::intersect(srng, sprng))
  uniqueLengthQuery <- sum(setDiffQ)
  uniqueLengthSubject <- sum(interesectS)
  
  queryElementsOutsideMaxDist <- sum(setDiffQ >= maxDist)
  subjectElementsOutsideMaxDist <- sum(interesectS >= maxDist)
  compatible <- (queryElementsOutsideMaxDist == 0) &
    (subjectElementsOutsideMaxDist == 0)
  DataFrame(uniqueLengthQuery, uniqueLengthSubject, compatible,
            queryElementsOutsideMaxDist, subjectElementsOutsideMaxDist)
}



#' The following function is implemented in R (GenomicAlignments), 
#' I just included the "within" option to make them significantly
#' faster+memorey friendly for this purpose (original code modified 
#' from the GenomicAlignments package, Author: Hervé Pagès, Valerie Obenchain,
#'  Martin Morgan)
# License Artistic-2.0
# https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments

#' @noRd
findSpliceOverlapsQuick <- function(query, subject, ignore.strand = FALSE) {
  olap <- findOverlaps(query, subject, ignore.strand = ignore.strand,
                       type = "within")
  olapEqual <- findOverlaps(query, subject, ignore.strand = ignore.strand,
                            type = "equal")
  if (length(olap) == 0L)
    return(GenomicAlignments:::.result(olap))
  
  query <- query[queryHits(olap)]
  subject <- subject[subjectHits(olap)]
  splice <- myGaps(query)
  
  compatible <- myCompatibleTranscription(query, subject, splice)
  strandSpecific <- all(strand(query) != "*")
  equal <- (!is.na(S4Vectors::match(olap, olapEqual)))
  unique <- myOneMatch(compatible, queryHits(olap))
  
  mcols(olap) <- DataFrame(compatible, equal, unique, strandSpecific)
  return(olap)
}

#' @param query query
#' @param subject subject
#' @param splice splice
#' @noRd
myCompatibleTranscription <- function(query, subject, splice) {
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
myOneMatch <- function(x, idx) {
  # License note: This function is adopted from the GenomicAlignments
  # package (Author: Hervé Pagès, Valerie Obenchain, Martin Morgan)
  # https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments
  xcnt <- rowsum(as.integer(x), idx)[, 1]
  oneMatch <- rep((xcnt == 1L), table(idx))
  unname(x & oneMatch)
}

#' @param motif motif
#' @noRd
spliceStrand <- function(motif) {
  NATURAL_INTRON_MOTIFS_RC <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAStringSet(GenomicAlignments::NATURAL_INTRON_MOTIFS)))
  
  motifStrand <- ifelse(motif %in% GenomicAlignments::NATURAL_INTRON_MOTIFS,
                        "+", "*")
  motifStrand[motif %in% NATURAL_INTRON_MOTIFS_RC] <- "-"
  return(motifStrand)
}


#' Function to reduce the start end end of the first and last elements in a 
#' granges list objects to a single basepair, helper to identify overlaps
#' based on splicing only (allow for flexible TSS/TES)
#' @param grangesList grangesList
#' @noRd
cutStartEndFromGrangesList <- function(grangesList) {
  unlistedExons <- unlist(grangesList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(grangesList)),
                                    names = NULL)
  startExonsSet <- (which((unlistedExons$exon_rank == 1 &
                             as.character(strand(unlistedExons)) != "-") | 
                            (unlistedExons$exon_endRank == 1 & 
                               as.character(strand(unlistedExons)) == "-")))
  endExonsSet <- (which((unlistedExons$exon_rank == 1 & 
                           as.character(strand(unlistedExons)) == "-") | 
                          (unlistedExons$exon_endRank == 1 
                           & as.character(strand(unlistedExons)) != "-")))
  
  start(unlistedExons[startExonsSet]) <- end(unlistedExons[startExonsSet]) - 1
  end(unlistedExons[endExonsSet]) <- start(unlistedExons[endExonsSet]) + 1
  
  return(relist(unlistedExons, partitioning))
}

#' @param grangesList grangesList
#' @param by defaults to 5
#' @noRd
extendGrangesListElements <- function(grangesList, by = 5) {
  unlistedExons <- unlist(grangesList, use.names = FALSE)
  partitioning <-
    PartitioningByEnd(cumsum(elementNROWS(grangesList)), names = NULL)
  start(unlistedExons) <- pmax(1, start(unlistedExons) - by)
  end(unlistedExons) <-
    pmin(seqlengths(unlistedExons)[as.character(seqnames(unlistedExons))],
         end(unlistedExons) + by, na.rm = TRUE)
  return(relist(unlistedExons, partitioning))
}


#' @param grangesList grangesList
#' @param minWidth defaults to 5
#' @param cutStartEnd defaults to FALSE
#' @noRd
dropGrangesListElementsByWidth <- function(grangesList, minWidth = 5,
                                           cutStartEnd = FALSE) {
  unlistedExons <- unlist(grangesList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(sum(width(grangesList) >= 
                                                 minWidth)), names = NULL)
  exonWidth <- width(unlistedExons)
  if (cutStartEnd) {
    startExonsSet <- (which((unlistedExons$exon_rank == 1 &
                               as.character(strand(unlistedExons)) != "-") |
                              (unlistedExons$exon_endRank == 1 &
                                 as.character(strand(unlistedExons)) == "-")))
    endExonsSet <- (which((unlistedExons$exon_rank == 1 &
                             as.character(strand(unlistedExons)) == "-") |
                            (unlistedExons$exon_endRank == 1 &
                               as.character(strand(unlistedExons)) != "-")))
    
    start(unlistedExons[startExonsSet]) <-
      end(unlistedExons[startExonsSet]) - 1
    end(unlistedExons[endExonsSet]) <-
      start(unlistedExons[endExonsSet]) + 1
  }
  unlistedExons <- unlistedExons[exonWidth >= minWidth]
  return(relist(unlistedExons, partitioning))
}

#' Function that selects the first exon from an IRangesList object
#' @param range IRangesList with elements required to be ordered by coordinates
#' @param stand strand
#' @noRd
selectStartExonFromRangesList <- function(range, strand){
  lastExons <- as.numeric(cumsum(elementNROWS(range)))
  startExonsSet <- c(1, lastExons[-(length(lastExons))]+1)
  startExonsSet[strand == "-"] <- lastExons[strand == "-"]
  return(unlist(range, use.names = FALSE)[startExonsSet])
}

#' Function that selects the last exon from an IRangesList object
#' @param range IRangesList with elements required to be ordered by coordinates
#' @param stand strand
#' @noRd
selectEndExonFromRangesList <- function(range, strand){
  endExonsSet <- as.numeric(cumsum(elementNROWS(range)))
  firstExons <- c(1, endExonsSet[-(length(endExonsSet))]+1)
  endExonsSet[strand == "-"] <- firstExons[strand == "-"]
  return(unlist(range, use.names = FALSE)[endExonsSet])
}

#' Function that selects the first N exons from a grangeslist object
#' (exon_rank is required)
#' @param grangesList grangesList
#' @param exonNumber defaults to 2
#' @noRd
selectStartExonsFromGrangesList <- function(grangesList, exonNumber = 2) {
  unlisted_granges <- unlist(grangesList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(pmin(elementNROWS(grangesList),
                                                exonNumber)), names = NULL)
  startExonsSet <- which(unlisted_granges$exon_rank <= exonNumber)
  return(relist(unlisted_granges[startExonsSet], partitioning))
}

#' Function that selects the last N exons from a grangeslist object 
#' (exon_endRank is required)
#' @describeIn selectStartExonsFromGrangesList grangesList
#' @noRd
selectEndExonsFromGrangesList <- function(grangesList, exonNumber = 2) {
  unlisted_granges <- unlist(grangesList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(pmin(elementNROWS(grangesList),
                                                exonNumber)), names = NULL)
  endExonsSet <- which(unlisted_granges$exon_endRank <= exonNumber)
  return(relist(unlisted_granges[endExonsSet], partitioning))
}
