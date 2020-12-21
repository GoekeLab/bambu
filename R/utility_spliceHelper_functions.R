# Note: several of the functions in this file are adopted from 
# the GenomicAlignments package (Author: Hervé Pagès, Valerie Obenchain,
# Martin Morgan)
# License Artistic-2.0
# https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments

## examples for test purposes
# query=rowRanges(seBambu.core)[c(
# 'ENST00000344579', # exon skipping, alternative TSS (-48), +, ENSG00000158109
# 'ENST00000270792', #intron retention subject 1(last exon),alt.TSS,alt.TES, +,
# 'ENST00000410032', # alternative first exon, exon skipping query: 2, 
# #exon skipping subject: 0, alternative TSS (2bp only), internalFirstExon.subject +
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
#     )]
  # subject=rowRanges(seBambu.core)[c('ENST00000378344',
  #                                   'ENST00000319041',
  #                                   'ENST00000338530',
  #                                   'ENST00000338530',
  #                                   'ENST00000338530',
  #                                   'ENST00000532718',
  #                                   'ENST00000263331',
  #                                   'ENST00000602866',
  #                                   'ENST00000585361')]
  # query <- rep(query,2000)
  # subject <- rep(subject,2000)
## To be modified after merging Yuk Kei's PR
annotateSpliceOverlapsByDist <- function(query, subject) {
  queryStartRng <- ranges(selectStartExonsFromGrangesList(query, exonNumber = 1))
  subjectStartRng <- ranges(selectStartExonsFromGrangesList(subject, exonNumber = 1))
  queryEndRng <- ranges(selectEndExonsFromGrangesList(query, exonNumber = 1))
  subjectEndRng <- ranges(selectEndExonsFromGrangesList(subject, exonNumber = 1))
  subjectFullRng <- ranges(subject)
  queryFullRng <- ranges(query)
  qSpRng <- ranges(myGaps(query))
  qSpRng[elementNROWS(qSpRng)==0] <- IRanges(start=1,end=1) # add mock intron
  sSpRng <- ranges(myGaps(subject))
  sSpRng[elementNROWS(sSpRng)==0] <- IRanges(start=1,end=1)# add mock intron

  ## start end 
  startEndTable <- tibble(queryId = names(query),
                          subjectId = names(subject),
                          start.first.query = start(unlist(queryStartRng)),
                          end.first.query = end(unlist(queryStartRng)),
                          start.last.query = start(unlist(queryEndRng)),
                          end.last.query = end(unlist(queryEndRng)),
                          start.first.subject = start(unlist(subjectStartRng)),
                          end.first.subject = end(unlist(subjectStartRng)),
                          start.last.subject = start(unlist(subjectEndRng)),
                          end.last.subject = end(unlist(subjectEndRng)),
                          strand = as.character(unlist(unique(strand(query)))))

  # calculate alternative First/last exons
  startEndTable <- startEndTable %>% mutate(alternativeFirstExon=!ifelse(strand!='-',
                                                                          start.first.query <= end.first.subject & end.first.query >= start.first.subject,
                                                                          end.first.query >= start.first.subject & start.first.query <= end.first.subject),
                                             alternativeLastExon=!ifelse(strand!='-',
                                                                         end.last.query >= start.last.subject & start.last.query <= end.last.subject,
                                                                         start.last.query <= end.last.subject & end.last.query >= start.last.subject),
                                             alternativeTSS=ifelse(strand!='-',start.first.subject-start.first.query,
                                                                  end.first.query- end.first.subject )*!alternativeFirstExon,
                                             alternativeTES=ifelse(strand!='-', 
                                                                   end.last.query-end.last.subject,
                                                                   start.last.subject-start.last.query)*!alternativeLastExon)

  startEndTable <- startEndTable %>% select(queryId,
                                            subjectId, 
                                            strand,
                                            alternativeFirstExon,
                                            alternativeLastExon,
                                            alternativeTSS,
                                            alternativeTES)


  ## annotate internal start and end first exons ##
  subjectList <- unlist(subjectFullRng)
  queryList <- unlist(queryFullRng)
  sSpRngList <- unlist(sSpRng)
  qSpRngList <- unlist(qSpRng)

  queryStart.subjectFull.Rng <- rep(unlist(queryStartRng),
                                    elementNROWS(subjectFullRng))
  queryStart.subjectFull.IdMap <- rep(1:length(queryStartRng),
                                      elementNROWS(subjectFullRng))

  subjectStart.queryFull.Rng <- rep(unlist(subjectStartRng),
                                    elementNROWS(queryFullRng))
  subjectStart.queryFull.IdMap <- rep(1:length(subjectStartRng),
                                      elementNROWS(queryFullRng))

  queryEnd.subjectFull.Rng <- rep(unlist(queryEndRng),
                                  elementNROWS(subjectFullRng))
  queryEnd.subjectFull.IdMap <- rep(1:length(queryEndRng),
                                    elementNROWS(subjectFullRng))

  subjectEnd.queryFull.Rng <- rep(unlist(subjectEndRng),
                                  elementNROWS(queryFullRng))
  subjectEnd.queryFull.IdMap <- rep(1:length(subjectEndRng),
                                    elementNROWS(queryFullRng))

  #for intron retention/exon skipping
  querySplice.subjectFullQuerySplice.Rng <- rep(unlist(qSpRng, use.names=F),
                                     rep(elementNROWS(subjectFullRng),
                                         times=elementNROWS(qSpRng)))
  subjectFullQuerySplice.Rng <- unlist(rep(subjectFullRng,
                                                       times= elementNROWS(qSpRng)), use.names=F)
  querySplice.subjectFullQuerySplice.IdMap <- rep(1:length(qSpRng),elementNROWS(qSpRng)*elementNROWS(subjectFullRng))

  subjectSplice.queryFullsubjectSplice.Rng <- rep(unlist(sSpRng, use.names=F),rep(elementNROWS(queryFullRng), times=elementNROWS(sSpRng)))
  queryFullSubjectSplice.Rng <- unlist(rep(queryFullRng, times= elementNROWS(sSpRng)), use.names=F)
  subjectSplice.queryFullSubjectSplice.IdMap <- rep(1:length(sSpRng), elementNROWS(sSpRng)*elementNROWS(queryFullRng))

  #for exon skipping

  subjectStart.querySplice.Rng <- rep(unlist(subjectStartRng),elementNROWS(qSpRng))
  subjectEnd.querySplice.Rng <- rep(unlist(subjectEndRng),elementNROWS(qSpRng))
  subjectStartEnd.querySplice.IdMap <- rep(1:length(subjectStartRng),elementNROWS(qSpRng))


  queryStart.subjectSplice.Rng <- rep(unlist(queryStartRng),elementNROWS(sSpRng))
  queryEnd.subjectSplice.Rng <- rep(unlist(queryEndRng),elementNROWS(sSpRng))
 queryStartEnd.subjectSplice.IdMap <- rep(1:length(queryStartRng),elementNROWS(sSpRng))




  #internal start query
  queryExonIntersect <- pintersect(queryStart.subjectFull.Rng,subjectList, resolve.empty='start.x')
  startEndTable$internalFirstExon.query <- tapply(width(queryExonIntersect), queryStart.subjectFull.IdMap, sum)!= 0 & startEndTable$alternativeFirstExon
  #internal start subject
  subjectExonIntersect <- pintersect(subjectStart.queryFull.Rng,queryList, resolve.empty='start.x')
  startEndTable$internalFirstExon.subject <- tapply(width(subjectExonIntersect), subjectStart.queryFull.IdMap, sum)!= 0 & startEndTable$alternativeFirstExon
  #internal end query
  queryExonIntersect=pintersect(queryEnd.subjectFull.Rng,subjectList, resolve.empty='start.x')
  startEndTable$internalLastExon.query <- tapply(width(queryExonIntersect), queryEnd.subjectFull.IdMap, sum)!=0 & startEndTable$alternativeLastExon
  #internal end subject
  subjectExonIntersect=pintersect(subjectEnd.queryFull.Rng,queryList, resolve.empty='start.x')
  startEndTable$internalLastExon.subject <- tapply(width(subjectExonIntersect), subjectEnd.queryFull.IdMap, sum)!=0 & startEndTable$alternativeLastExon
  #intron retention subject
  subjectIntronRetention <-punion(querySplice.subjectFullQuerySplice.Rng, subjectFullQuerySplice.Rng, fill.gap=TRUE) == subjectFullQuerySplice.Rng
  startEndTable$intronRetention.subject <- tapply(subjectIntronRetention, querySplice.subjectFullQuerySplice.IdMap, sum)
  #intron retention query
  queryIntronRetention <-punion(subjectSplice.queryFullsubjectSplice.Rng, queryFullSubjectSplice.Rng, fill.gap=TRUE) == queryFullSubjectSplice.Rng
  startEndTable$intronRetention.query <- tapply(queryIntronRetention, subjectSplice.queryFullSubjectSplice.IdMap, sum)


  #exon skipping query
  queryExonSkipping <- punion(querySplice.subjectFullQuerySplice.Rng,
                              subjectFullQuerySplice.Rng,
                              fill.gap=TRUE) == querySplice.subjectFullQuerySplice.Rng 

  queryFirstExonInIntron <-punion(qSpRngList, 
                                  subjectStart.querySplice.Rng, 
                                  fill.gap=TRUE) == qSpRngList

  queryLastExonInIntron <-punion(qSpRngList, 
                                 subjectEnd.querySplice.Rng, 
                                 fill.gap=TRUE) == qSpRngList

  startEndTable$exonSkipping.query <- pmax(0,
                                           tapply(queryExonSkipping, 
                                                  querySplice.subjectFullQuerySplice.IdMap, 
                                                  sum) -  
                                             tapply(queryFirstExonInIntron, 
                                                    subjectStartEnd.querySplice.IdMap, 
                                                    sum) - 
                                             tapply(queryLastExonInIntron, 
                                                    subjectStartEnd.querySplice.IdMap, 
                                                    sum))

  #exon skipping subject

  subjectExonSkipping <- punion(subjectSplice.queryFullsubjectSplice.Rng,
                                queryFullSubjectSplice.Rng, 
                                fill.gap=TRUE) == subjectSplice.queryFullsubjectSplice.Rng
  subjectFirstExonInIntron <- punion(sSpRngList, 
                                     queryStart.subjectSplice.Rng, 
                                     fill.gap=TRUE) == sSpRngList
  subjectLastExonInIntron <- punion(sSpRngList, 
                                    queryEnd.subjectSplice.Rng, 
                                    fill.gap=TRUE) == sSpRngList


  startEndTable$exonSkipping.subject <- pmax(0, 
                                             tapply(subjectExonSkipping,
                                                    subjectSplice.queryFullSubjectSplice.IdMap,
                                                    sum) -
                                               tapply(subjectFirstExonInIntron, 
                                                      queryStartEnd.subjectSplice.IdMap,
                                                      sum) - 
                                               tapply(subjectLastExonInIntron, 
                                                      queryStartEnd.subjectSplice.IdMap,
                                                      sum))

  # exon 3' splice site

  queryExonEnd <- start(querySplice.subjectFullQuerySplice.Rng)<end(subjectFullQuerySplice.Rng) & 
    start(querySplice.subjectFullQuerySplice.Rng)>start(subjectFullQuerySplice.Rng) &
    end(querySplice.subjectFullQuerySplice.Rng)>end(subjectFullQuerySplice.Rng)
  queryExonStart <- end(querySplice.subjectFullQuerySplice.Rng)<end(subjectFullQuerySplice.Rng) & 
    end(querySplice.subjectFullQuerySplice.Rng)>start(subjectFullQuerySplice.Rng) &
    start(querySplice.subjectFullQuerySplice.Rng)<start(subjectFullQuerySplice.Rng)

  exon5Prime <- tapply(queryExonStart, querySplice.subjectFullQuerySplice.IdMap, sum)
  exon3Prime <- tapply(queryExonEnd, querySplice.subjectFullQuerySplice.IdMap, sum)

  subjectStartExonStartExtension <- start(subjectStart.querySplice.Rng)<end(qSpRngList) & end(subjectStart.querySplice.Rng)>end(qSpRngList) & start(qSpRngList)< start(subjectStart.querySplice.Rng)
  subjectStartExonStartExtension <- tapply(subjectStartExonStartExtension, subjectStartEnd.querySplice.IdMap, sum)

  subjectStartExonEndExtension <- end(subjectStart.querySplice.Rng)>start(qSpRngList) & start(subjectStart.querySplice.Rng)<start(qSpRngList) & end(qSpRngList)> end(subjectStart.querySplice.Rng)
  subjectStartExonEndExtension <- tapply(subjectStartExonEndExtension, subjectStartEnd.querySplice.IdMap, sum)

  subjectEndExonStartExtension <- start(subjectEnd.querySplice.Rng)<end(qSpRngList) & end(subjectEnd.querySplice.Rng)>end(qSpRngList) & start(qSpRngList)< start(subjectEnd.querySplice.Rng)
  subjectEndExonStartExtension <- tapply(subjectEndExonStartExtension, subjectStartEnd.querySplice.IdMap, sum)

  subjectEndExonEndExtension <- end(subjectEnd.querySplice.Rng)>start(qSpRngList) & start(subjectEnd.querySplice.Rng)<start(qSpRngList) & end(qSpRngList)> end(subjectEnd.querySplice.Rng)
  subjectEndExonEndExtension <- tapply(subjectEndExonEndExtension, subjectStartEnd.querySplice.IdMap, sum)

  startEndTable$exon5Prime <- exon5Prime
  startEndTable$exon5Prime[startEndTable$strand!='-'] <- startEndTable$exon5Prime[startEndTable$strand!='-'] - subjectStartExonStartExtension[startEndTable$strand!='-']
  startEndTable$exon5Prime[startEndTable$strand=='-'] <- exon3Prime[startEndTable$strand=='-'] - subjectStartExonEndExtension[startEndTable$strand=='-']

  startEndTable$exon3Prime <- exon3Prime
  startEndTable$exon3Prime[startEndTable$strand!='-'] <- startEndTable$exon3Prime[startEndTable$strand!='-'] - subjectEndExonEndExtension[startEndTable$strand!='-']
  startEndTable$exon3Prime[startEndTable$strand=='-'] <- exon5Prime[startEndTable$strand=='-'] - subjectEndExonStartExtension[startEndTable$strand=='-']
  return(startEndTable)
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
