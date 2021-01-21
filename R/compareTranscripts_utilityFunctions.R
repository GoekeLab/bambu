# This file contains functions that are based on ranges (rather than GRanges)
# with the exception of getStrandFromGrList

#' extract strand from GRangesList
#' @description This function takes a GRangesList and
#' returns a vector with the strand for each list entry.
#' This function assumes that all elements for each list
#' entry have the same strand.
#' @usage getStrangeFromGrList(grl)
#' @param grl a GRangesList
#' @return an Rle object with strand information
#' @examples
#' query <- readRDS(system.file("extdata", 
#'     "annotateSpliceOverlapByDist_testQuery.rds",
#'     package = "bambu"))
#' strand <- as.character(getStrandFromGrList(query))
#' @noRd
getStrandFromGrList <- function(grl) { 
    return(unlist(strand(grl), use.names = FALSE)[cumsum(elementNROWS(grl))]) 
}


#' Function that selects the first/last exon from an IRangesList object
#' @param range IRangesList with elements required to be ordered by coordinates
#' @param stand strand
#' @noRd
selectStartEndExonFromRangesList <- function(range, strand, direction){
    exons <- as.numeric(cumsum(elementNROWS(range)))
    exonsSet <- c(1, exons[-(length(exons))] + 1)
    if (direction == "start") {
        exonsSet[strand == "-"] <- exons[strand == "-"]
        return(unlist(range, use.names = FALSE)[exonsSet])
    } else{
        exons[strand == "-"] <- exonsSet[strand == "-"] 
        return(unlist(range, use.names = FALSE)[exons])
    }
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
        (get(direction)(subjectTerminalExonRng) - 
        get(direction)(queryTerminalExonRng))
    alternativeTerminal[strand == "-"] <- (-1)^(direction == "end") * 
        (get(setdiff(direction_names, direction))(
            queryTerminalExonRng[strand == '-']) -
        get(setdiff(direction_names, direction))(
            subjectTerminalExonRng[strand == '-']))
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
        mcols(exon.Full.Rng)$matchRng, resolve.empty = 'start.x')
    internalStartEndVector <- tapply(width(exonIntersect),
        mcols(exon.Full.Rng)$IdMap, sum) != 0 & alternativeFirstLastExon
    return(internalStartEndVector)
}

#' annotate intron retention
#' @description This function checks whether
#' there is intron retention by overlapping
#' the intron ranges of matching transcripts.
#' @noRd
annotateIntronRetent <- function(spliceRng, fullRng){
    splice.FullSplice.Rng <- expandRangesList(spliceRng, fullRng)
    intronRetention <- punion(splice.FullSplice.Rng, 
        mcols(splice.FullSplice.Rng)$matchRng, fill.gap = TRUE) ==
    mcols(splice.FullSplice.Rng)$matchRng
    intronRetentionVector <- tapply(intronRetention,
        mcols(splice.FullSplice.Rng)$IdMap, sum)
    return(intronRetentionVector)
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
        fill.gap = TRUE) == splice.FullSplice.Rng
    firstExonInIntron <- punion(mcols(start.Splice.Rng)$matchRng, 
        start.Splice.Rng, fill.gap = TRUE) == mcols(start.Splice.Rng)$matchRng
    lastExonInIntron <- punion(mcols(start.Splice.Rng)$matchRng, end.Splice.Rng,
        fill.gap = TRUE) == mcols(start.Splice.Rng)$matchRng
    exonSkippingVector <- pmax(0, tapply(exonSkipping,
        mcols(splice.FullSplice.Rng)$IdMap, sum) -
        tapply(firstExonInIntron, 
        mcols(start.Splice.Rng)$IdMap, sum) - 
        tapply(lastExonInIntron, 
        mcols(start.Splice.Rng)$IdMap, sum))
    return(exonSkippingVector)
}


#' annotate exon splicing
#' @description This function checks whether 
#' there is alternative splicing in the 5'/3'
#' end of an exon.
#' @importFrom dplyr tibble %>% select
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
    exonEnd <- startMatch & !endMatch
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
    exStrandNeg <- exonSplicingTable$strand == '-'
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
    exonSplicingTable <- exonSplicingTable %>% 
        dplyr::select(exon5Prime,exon3Prime)
    return(exonSplicingTable)
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
    exonStartExtension <- splice.start < match.startSplice.end &
        splice.end > match.startSplice.end &
        match.startSplice.start < splice.start
    exonStartExtension <- tapply(exonStartExtension, spliceIdMap, sum)
    return(exonStartExtension)
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
    exonEndExtension <- splice.end > match.startSplice.start &
        splice.start < match.startSplice.start &
        match.startSplice.end > splice.end
    exonEndExtension <- tapply(exonEndExtension, spliceIdMap, sum)
    return(exonEndExtension)
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
    mcols(processedRng)$IdMap <- rep(seq_along(ranges),elementNROWS(target))
    mcols(processedRng)$matchRng <- unlist(target, use.names = FALSE)
    return(processedRng)
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
    processedRng <- rep(unlist(rglist, use.names = FALSE),
        rep(elementNROWS(target),times = elementNROWS(rglist)))
    mcols(processedRng)$IdMap <- rep(seq_along(rglist),
        elementNROWS(rglist) * elementNROWS(target))
    mcols(processedRng)$matchRng <- unlist(rep(target, 
        times = elementNROWS(rglist)), use.names = FALSE)
    return(processedRng)
}
