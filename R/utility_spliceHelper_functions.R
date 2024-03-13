# Note: several of the functions in this file are adopted from 
# the GenomicAlignments package (Author: Hervé Pagès, Valerie Obenchain,
# Martin Morgan)
# License Artistic-2.0
# https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments


#' calculate distance between first and last exon matches
#' @param queryExon a query start or end exon ranges
#' @param subjectExon a subject start or end exon ranges
#' @param subjectFull a full subject ranges object
#' @param subjectList a full subject list
#' @importFrom GenomicRanges setdiff
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
#' @importFrom S4Vectors match
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
    subject = split_olaps(subject, olap)
    #subject <- subject[subjectHits(olap)]
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

split_olaps = function(grl,ov, ncore = 16){
    subjects = parallel::mclapply(split(subjectHits(ov),cut(1:length(subjectHits(ov)),100)),function(subjectHit){
        grl[subjectHit]
    },mc.cores = ncore)
    subjects = subjects[lengths(subjects)!=0]
    subjects = do.call("c",unlist(subjects, use.names = FALSE))  
}

#' check whether error with start sequence
#' @noRd
checkStartSequence <- function(olap, firstLastSeparate, queryStart,
        subjectStart, queryEnd,subjectEnd, subjectFull, subjectList){
    if (length(olap)) {
        qHits <- queryHits(olap)
        subHits <- subjectHits(olap)
        queryStart <- ranges(queryStart[qHits])
        subjectStart <- ranges(subjectStart[subHits])
        queryEnd <- ranges(queryEnd[qHits])
        subjectEnd <- ranges(subjectEnd[subHits])
        subjectFull <- ranges(subjectFull[subHits])
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
#' @importFrom methods as
#' @importFrom GenomicRanges GRanges setdiff 
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
        insert_gaps <- as(ranges(.insertGaps(x)), "CompressedIRangesList")
        gaps <- setdiff(gaps, insert_gaps)
    }
    
    idx <- elementNROWS(gaps) != 0
    ## FIXME : can't handle lists with empty elements
    ##         'start' and 'end' not quite right here
    firstseg <- start(PartitioningByWidth(x))
    seqnms <- rep(seqnames(flat)[firstseg], elementNROWS(gaps))
    strand <- rep(strand(flat)[firstseg], elementNROWS(gaps))
    gr <- relist(GRanges(seqnms, unlist(gaps, use.names = FALSE), strand), gaps)
    return(gr)
    } else {
    ### FIXME: does not handle query.break column yet
    return(GenomicRanges::setdiff(range(x), x))
    }
}
# myGaps <- .GenomicAlignments:::.gaps

#' @noRd
.isNumericOrNAs <- S4Vectors:::isNumericOrNAs


#' @importFrom GenomicRanges setdiff intersect 
#' @noRd
rangesDist <- function(query, subject, splice, maxDist) {
    qrng <- ranges(query)
    srng <- ranges(subject)
    sprng <- ranges(splice)
    
    setDiffQ <- width(setdiff(qrng, srng))
    interesectS <- width(intersect(srng, sprng))
    uniqueLengthQuery <- sum(setDiffQ)
    uniqueLengthSubject <- sum(interesectS)
    
    queryElementsOutsideMaxDist <- sum(setDiffQ >= maxDist)
    subjectElementsOutsideMaxDist <- sum(interesectS >= maxDist)
    compatible <- (queryElementsOutsideMaxDist == 0) &
        (subjectElementsOutsideMaxDist == 0)
    return(DataFrame(uniqueLengthQuery, uniqueLengthSubject, compatible,
        queryElementsOutsideMaxDist, subjectElementsOutsideMaxDist))
}



#' The following function is implemented in R (GenomicAlignments), 
#' I just included the "within" option to make them significantly
#' faster+memorey friendly for this purpose (original code modified 
#' from the GenomicAlignments package, Author: Hervé Pagès, Valerie Obenchain,
#'  Martin Morgan)
# License Artistic-2.0
# https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments

#' @importFrom S4Vectors match
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
#' @importFrom GenomicRanges setdiff intersect 
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



#' Function to reduce the start end end of the first and last elements in a 
#' granges list objects to a single basepair, helper to identify overlaps
#' based on splicing only (allow for flexible TSS/TES)
#' @param grangesList grangesList
#' @noRd
cutStartEndFromGrangesList <- function(grangesList) {
    unlistedExons <- unlist(grangesList, use.names = FALSE)
    partitioning <- PartitioningByEnd(cumsum(elementNROWS(grangesList)),
                                    names = NULL)
    exonRank <- unlistedExons$exon_rank
    exonEndRank <- unlistedExons$exon_endRank
    exonStrand <- as.character(strand(unlistedExons))
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
