
#' Extend annotations
#' @inheritParams bambu
#' @noRd
isore.extendAnnotations <- function(se, annotationGrangesList,
                                    remove.subsetTx = TRUE, min.readCount = 2, 
                                    min.readFractionByGene = 0.05, min.sampleNumber = 1, min.exonDistance = 35, 
                                    min.exonOverlap = 10, min.primarySecondaryDist = 5,
                                    min.primarySecondaryDistStartEnd = 5, prefix = "", verbose = FALSE){
    start.ptm <- proc.time() # timer for verbose output
    filterSet1 <- FALSE
    if (nrow(se) > 0) filterSet1 <-
        rowSums(assays(se)$counts >= min.readCount) >= min.sampleNumber
    if (sum(filterSet1) > 0) {
        if (any(rowData(se)$confidenceType == "highConfidenceJunctionReads" &
                filterSet1)) { ## (1) Spliced Reads
            byReadClass <- exonsintronsByReadClass(
                se, annotationGrangesList, filterSet1)
            seFilteredSpliced <- extdannotateSplicedReads(
                byReadClass$exons, byReadClass$intron, annotationGrangesList,
                byReadClass$seFilteredSpliced, min.exonDistance, verbose,
                min.primarySecondaryDist, min.primarySecondaryDistStartEnd)
            end.ptm <- proc.time()
            if (verbose)
                message("extended annotations for spliced reads in ",
                        round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
        }
        ## unspliced transcripts
        SEnRng <- extdannotateUnsplicedReads(
            se, seFilteredSpliced, byReadClass$exons, annotationGrangesList,
            filterSet1, min.exonOverlap, verbose)
        
        seCombined <- SEnRng$seCombined
        exonRangesCombined <- SEnRng$exonRangesCombined
        # assign gene IDs based on exon match
        seCombined <- assignGeneIDexonMatch(
            seCombined, exonRangesCombined, annotationGrangesList,
            min.exonOverlap, prefix, verbose)
        ## filter out transcripts
        extendedAnnotationRanges <- filterTranscripts(
            seCombined, annotationGrangesList, exonRangesCombined, prefix,
            min.readFractionByGene, min.sampleNumber, remove.subsetTx, verbose)
        return(extendedAnnotationRanges)
    } else {
        return(annotationGrangesList)
    }
}


#' generate exon/intron ByReadClass objects
#' @noRd
exonsintronsByReadClass <- function(se, annotationGrangesList, filterSet1){
    seFilteredSpliced <- se[rowData(se)$confidenceType ==
                                "highConfidenceJunctionReads" & filterSet1, ]
    mcols(seFilteredSpliced)$GENEID <- NA
    intronsByReadClass <- GenomicRanges::makeGRangesListFromFeatureFragments(
        seqnames = rowData(seFilteredSpliced)$chr,
        fragmentStarts = rowData(seFilteredSpliced)$intronStarts,
        fragmentEnds = rowData(seFilteredSpliced)$intronEnds,
        strand = rowData(seFilteredSpliced)$strand)
    names(intronsByReadClass) <- seq_along(intronsByReadClass)
    GenomeInfoDb::seqlevels(intronsByReadClass) <-
        unique(c(GenomeInfoDb::seqlevels(intronsByReadClass), 
                 GenomeInfoDb::seqlevels(annotationGrangesList)))
    exonsByReadClass <- createExonByReadClass(
        seFilteredSpliced, annotationGrangesList
    )
    return(list("exons" = exonsByReadClass, 
                "introns" = intronsByReadClass,
                "seFilteredSpliced" = seFilteredSpliced))
}


#' create exonsByReadClass
#' @param seFilteredSpliced a SummarizedExperiment object
#' for filtered spliced reads
#' @param annotationGrangesList annotation GRangesList object
#' @noRd
createExonByReadClass <- function(seFilteredSpliced, annotationGrangesList) {
    exonEndsShifted <- paste(rowData(seFilteredSpliced)$intronStarts,
                             as.integer(rowData(seFilteredSpliced)$end + 1), sep = ",")
    exonStartsShifted <- paste(as.integer(rowData(seFilteredSpliced)$start - 1),
                               rowData(seFilteredSpliced)$intronEnds, sep = ",")
    exonsByReadClass <- GenomicRanges::makeGRangesListFromFeatureFragments(
        seqnames = rowData(seFilteredSpliced)$chr,
        fragmentStarts = exonStartsShifted,
        fragmentEnds = exonEndsShifted,
        strand = rowData(seFilteredSpliced)$strand)
    exonsByReadClass <- narrow(exonsByReadClass, start = 2, end = -2)
    # correct junction to exon differences in coordinates
    names(exonsByReadClass) <- seq_along(exonsByReadClass)
    # add exon start and exon end rank
    unlistData <- unlist(exonsByReadClass, use.names = FALSE)
    partitioning <- PartitioningByEnd(cumsum(elementNROWS(exonsByReadClass)),
                                      names = NULL)
    exon_rank <- lapply(width((partitioning)), seq, from = 1)
    # * assumes positive for exon ranking
    exon_rank[which(rowData(seFilteredSpliced)$strand == "-")] <-
        lapply(exon_rank[which(rowData(seFilteredSpliced)$strand == "-")], rev) 
    exon_endRank <- lapply(exon_rank, rev)
    unlistData$exon_rank <- unlist(exon_rank)
    unlistData$exon_endRank <- unlist(exon_endRank)
    exonsByReadClass <- relist(unlistData, partitioning)
    GenomeInfoDb::seqlevels(exonsByReadClass) <-
        unique(c(GenomeInfoDb::seqlevels(exonsByReadClass), 
                 GenomeInfoDb::seqlevels(annotationGrangesList)))
    return(exonsByReadClass)
}


#' extended annotations for spliced reads
#' @param exonsByReadClass exonsByReadClass
#' @param intronsByReadClass intronsByReadClass
#' @param annotationGrangesList annotationGrangesList
#' @param seFilteredSpliced seFilteredSpliced
#' @param min.exonDistance min.exonDistance
#' @param verbose verbose
#' @noRd
extdannotateSplicedReads <- function(exonsByReadClass,intronsByReadClass, 
                                     annotationGrangesList, seFilteredSpliced, min.exonDistance, verbose,
                                     min.primarySecondaryDist, min.primarySecondaryDistStartEnd) {
    ovExon <- findSpliceOverlapsQuick(
        cutStartEndFromGrangesList(exonsByReadClass),
        cutStartEndFromGrangesList(annotationGrangesList))
    classificationTable <- 
        data.frame(matrix("", nrow = length(seFilteredSpliced), ncol = 9), 
                   stringsAsFactors = FALSE)
    colnames(classificationTable) <- c("equal", "compatible", "newWithin",
                                       "newLastJunction","newFirstJunction","newJunction", "allNew", 
                                       "newFirstExon", "newLastExon")
    classificationTable$equal[queryHits(ovExon[mcols(ovExon)$equal])[
        !duplicated(queryHits(ovExon[mcols(ovExon)$equal]))]] <- "equal"
    classificationTable$compatible[queryHits(ovExon[mcols(ovExon)$compatible])[
        !duplicated(queryHits(ovExon[mcols(ovExon)$compatible]))]] <- "compatible"
    classificationTable$compatible[classificationTable$equal == "equal"] <- ""
    ## annotate with transcript and gene Ids
    mcols(seFilteredSpliced)$GENEID[queryHits(ovExon[mcols(ovExon)$compatible][
        !duplicated(queryHits(ovExon[mcols(ovExon)$compatible]))])] <-
        mcols(annotationGrangesList[subjectHits(
            ovExon[mcols(ovExon)$compatible])[!duplicated(
                queryHits(ovExon[mcols(ovExon)$compatible]))]])$GENEID
    # annotate with compatible gene id,
    mcols(seFilteredSpliced)$GENEID[queryHits(ovExon[mcols(ovExon)$equal][
        !duplicated(queryHits(ovExon[mcols(ovExon)$equal]))])] <-
        mcols(annotationGrangesList[subjectHits(ovExon[mcols(ovExon)$equal])[
            !duplicated(queryHits(ovExon[mcols(ovExon)$equal]))]])$GENEID
    # annotate as identical, using intron matches
    unlistedIntrons <- unlist(intronsByReadClass, use.names = TRUE)
    partitioning <- PartitioningByEnd(cumsum(elementNROWS(intronsByReadClass)),
                                      names = NULL)
    unlistedIntronsAnnotations <- unlist(myGaps(annotationGrangesList))
    mcols(unlistedIntronsAnnotations)$GENEID <- mcols(annotationGrangesList
    )$GENEID[match(names(unlistedIntronsAnnotations),
                   mcols(annotationGrangesList)$TXNAME)]
    classificationTable <- 
        updateWIntronMatches(unlistedIntrons, unlistedIntronsAnnotations,
                             partitioning, classificationTable, annotationGrangesList,
                             seFilteredSpliced, exonsByReadClass, min.exonDistance,
                             min.primarySecondaryDist, min.primarySecondaryDistStartEnd)
    mcols(seFilteredSpliced)$readClassType <-
        apply(classificationTable, 1, paste, collapse = "")
    return(seFilteredSpliced)
}



#' update classificationTable
#' @noRd
updateWIntronMatches <- function(unlistedIntrons, unlistedIntronsAnnotations,
                                 partitioning, classificationTable, annotationGrangesList,
                                 seFilteredSpliced, exonsByReadClass, min.exonDistance,
                                 min.primarySecondaryDist, min.primarySecondaryDistStartEnd){
    intronMatches <- GenomicRanges::match(
        unlistedIntrons,unique(unlistedIntronsAnnotations),nomatch = 0) > 0
    intronMatchesList <- relist(intronMatches, partitioning)
    classificationTable$newWithin[all(intronMatchesList) & !(
        classificationTable$compatible == "compatible" |
            classificationTable$equal == "equal")] <- "newWithin"
    intronMatchListRev <- lapply(intronMatchesList, rev)
    lastJunctionMatch <- unlist(lapply(intronMatchListRev,`[[`, 1))
    firstJunctionMatch <- unlist(lapply(intronMatchesList, `[[`, 1))
    classificationTable$newLastJunction[which(rowData(
        seFilteredSpliced)$strand == "+" & !lastJunctionMatch &
            any(intronMatchesList) | (rowData(
                seFilteredSpliced)$strand == "-" & (!firstJunctionMatch &
                                                        any(intronMatchesList))))] <- "newLastJunction"
    classificationTable$newFirstJunction[which(rowData(
        seFilteredSpliced)$strand == "+" & !firstJunctionMatch &
            any(intronMatchesList) | (rowData(
                seFilteredSpliced)$strand == "-" & (!lastJunctionMatch &
                                                        any(intronMatchesList))))] <- "newFirstJunction"
    classificationTable$newJunction[(
        sum(!intronMatchesList) > !firstJunctionMatch + !lastJunctionMatch) &
            any(intronMatchesList)] <- "newJunction"
    classificationTable$allNew[!any(intronMatchesList)] <- "allNew"
    ## assign gene ids based on the max # of matching introns/splice junctions
    overlapsNewIntronsAnnotatedIntrons <- findOverlaps(unlistedIntrons,
                                                       unlistedIntronsAnnotations, type = "equal",
                                                       select = "all", ignore.strand = FALSE)
    if (length(overlapsNewIntronsAnnotatedIntrons)) {
        distNewTxByQuery <- assignGeneIDbyMaxMatch(
            unlistedIntrons,unlistedIntronsAnnotations,
            overlapsNewIntronsAnnotatedIntrons, exonsByReadClass,
            seFilteredSpliced, annotationGrangesList, min.exonDistance,
            min.primarySecondaryDist, min.primarySecondaryDistStartEnd)
        classificationTable$compatible[distNewTxByQuery$queryHits[
            distNewTxByQuery$compatible]] <- "compatible"
        classificationTable$newFirstExon[distNewTxByQuery$queryHits[
            !distNewTxByQuery$startMatch]] <- "newFirstExon"
        classificationTable$newFirstExon[
            classificationTable$newFirstJunction != "newFirstJunction"] <- ""
        classificationTable$newLastExon[distNewTxByQuery$queryHits[
            !distNewTxByQuery$endMatch]] <- "newLastExon"
        classificationTable$newLastExon[
            classificationTable$newLastJunction != "newLastJunction"] <- ""
    }
    return(classificationTable)
}


#' assign gene id by maximum match
#' @noRd
assignGeneIDbyMaxMatch <- function(unlistedIntrons,
                                   unlistedIntronsAnnotations, overlapsNewIntronsAnnotatedIntrons,
                                   exonsByReadClass, seFilteredSpliced, annotationGrangesList,
                                   min.exonDistance, min.primarySecondaryDist,
                                   min.primarySecondaryDistStartEnd) {
    maxGeneCountPerNewTx <- as_tibble(data.frame(txId =
                                                     names(unlistedIntrons)[queryHits(overlapsNewIntronsAnnotatedIntrons)],
                                                 geneId = mcols(unlistedIntronsAnnotations)$GENEID[
                                                     subjectHits(overlapsNewIntronsAnnotatedIntrons)],
                                                 stringsAsFactors = FALSE)) %>%
        group_by(txId, geneId) %>%
        summarise(geneCount = n()) %>%
        group_by(txId) %>%
        filter(geneCount == max(geneCount)) %>%
        filter(!duplicated(txId)) %>%
        ungroup()
    geneIdByIntron <- rep(NA, length(exonsByReadClass))
    geneIdByIntron <- maxGeneCountPerNewTx$geneId[
        match(names(exonsByReadClass), maxGeneCountPerNewTx$txId)]
    mcols(seFilteredSpliced)$GENEID[is.na(mcols(seFilteredSpliced)$GENEID)] <-
        geneIdByIntron[is.na(mcols(seFilteredSpliced)$GENEID)]
    distNewTx <- calculateDistToAnnotation(
        exByTx = exonsByReadClass,
        exByTxRef = annotationGrangesList,
        maxDist = min.exonDistance,
        primarySecondaryDist = min.primarySecondaryDist,
        primarySecondaryDistStartEnd = min.primarySecondaryDistStartEnd,
        ignore.strand = FALSE)
    distNewTxByQuery <- distNewTx %>%
        group_by(queryHits) %>%
        summarise(
            minDist = min(dist),
            startMatch = any(startMatch),
            endMatch = any(endMatch),
            compatible = any(compatible))
    ## note: here is more information that can be used to filter and annotate!
    return(distNewTxByQuery)
}


#' Calculate distance from read class to annotation
#' @param exByTx exByTx
#' @param exByTxRef exByTxRef
#' @param maxDist defaults to 35
#' @param primarySecondaryDist defaults to 5
#' @param ignore.strand defaults to FALSE
#' @noRd
calculateDistToAnnotation <- function(exByTx, exByTxRef, maxDist = 35,
                                      primarySecondaryDist = 5, primarySecondaryDistStartEnd = 5,
                                      ignore.strand = FALSE) {
    # (1)  find overlaps of read classes with annotated transcripts,
    spliceOverlaps <- findSpliceOverlapsByDist(exByTx, exByTxRef,
                                               maxDist = maxDist, firstLastSeparate = TRUE,
                                               dropRangesByMinLength = TRUE, cutStartEnd = TRUE,
                                               ignore.strand = ignore.strand)
    txToAnTableFiltered <- genFilteredAnTable(spliceOverlaps,
                                              primarySecondaryDist, DistCalculated = FALSE)
    # (2) calculate splice overlap for any not in the list (new exon >= 35bp)
    setTMP <- unique(txToAnTableFiltered$queryHits)
    spliceOverlaps_rest <- findSpliceOverlapsByDist(exByTx[-setTMP],
                                                    exByTxRef, maxDist = 0, type = "any", firstLastSeparate = TRUE,
                                                    dropRangesByMinLength = FALSE, cutStartEnd = TRUE,
                                                    ignore.strand = ignore.strand)
    txToAnTableRest <-
        genFilteredAnTable(spliceOverlaps_rest, primarySecondaryDist,
                           exByTx = exByTx, setTMP = setTMP, DistCalculated = FALSE)
    # (3) find overlaps for remaining reads 
    setTMPRest <- unique(c(txToAnTableRest$queryHits, setTMP))
    txToAnTableRestStartEnd <- NULL
    if (length(exByTx[-setTMPRest])) {
        spliceOverlaps_restStartEnd <-
            findSpliceOverlapsByDist(exByTx[-setTMPRest], exByTxRef,
                                     maxDist = 0, type = "any", firstLastSeparate = TRUE,
                                     dropRangesByMinLength = FALSE,
                                     cutStartEnd = FALSE, ignore.strand = ignore.strand)
        if (length(spliceOverlaps_restStartEnd)) {
            txToAnTableRestStartEnd <-
                genFilteredAnTable(spliceOverlaps_restStartEnd,
                                   primarySecondaryDist, exByTx = exByTx,
                                   setTMP = setTMPRest, DistCalculated = TRUE)
        }
    }
    txToAnTableFiltered <- rbind( txToAnTableFiltered,
                                  txToAnTableRest, txToAnTableRestStartEnd ) %>% ungroup()
    txToAnTableFiltered$readClassId <-
        names(exByTx)[txToAnTableFiltered$queryHits]
    txToAnTableFiltered$annotationTxId <-
        names(exByTxRef)[txToAnTableFiltered$subjectHits]
    return(txToAnTableFiltered)
}


#' generate filtered annotation table
#' @param spliceOverlaps an output from  findSpliceOverlapsByDist()
#' @param primarySecondaryDist default 5
#' @param primarySecondaryDistStartEnd default 5
#' @param exByTx default NULL
#' @param setTMP default NULL
#' @param DistCalculated default FALSE
#' @noRd
genFilteredAnTable <- function(spliceOverlaps, primarySecondaryDist = 5,
                               primarySecondaryDistStartEnd = 5, exByTx = NULL, setTMP = NULL,
                               DistCalculated = FALSE) {
    ## initiate the table
    if (isFALSE(DistCalculated)) {
        txToAnTable <- as_tibble(spliceOverlaps) %>% group_by(queryHits) %>%
            mutate(dist = uniqueLengthQuery + uniqueLengthSubject) %>%
            mutate(txNumber = n())
    } else {
        txToAnTable <- as_tibble(spliceOverlaps) %>% group_by(queryHits) %>%
            mutate(dist = uniqueLengthQuery + uniqueLengthSubject +
                       uniqueStartLengthQuery + uniqueEndLengthQuery) %>%
            mutate(txNumber = n())
    }
    ## change query hits for step 2 and 3
    if (!is.null(exByTx)) {
        txToAnTable$queryHits <-
            (seq_along(exByTx))[-setTMP][txToAnTable$queryHits]
    }
    ## todo: check filters, what happens to reads with only start and end match?
    if (isFALSE(DistCalculated)) {
        txToAnTableFiltered <- txToAnTable %>%
            group_by(queryHits) %>%
            arrange(queryHits, dist) %>%
            filter(dist <= (min(dist) + primarySecondaryDist)) %>%
            filter(queryElementsOutsideMaxDist + 
                       subjectElementsOutsideMaxDist == 
                       min(queryElementsOutsideMaxDist +
                               subjectElementsOutsideMaxDist)) %>% 
            filter((uniqueStartLengthQuery <= primarySecondaryDistStartEnd &
                        uniqueEndLengthQuery <= primarySecondaryDistStartEnd) ==
                       max(uniqueStartLengthQuery <=
                               primarySecondaryDistStartEnd & uniqueEndLengthQuery <=
                               primarySecondaryDistStartEnd)) %>%
            mutate(txNumberFiltered = n())
    } else {
        txToAnTableFiltered <- txToAnTable %>%
            group_by(queryHits) %>%
            arrange(queryHits, dist) %>%
            filter(dist <= (min(dist) + primarySecondaryDist)) %>%
            mutate(txNumberFiltered = n())
    }
    return(txToAnTableFiltered)
}

#' extended annotations for unspliced reads
#' @param se a summarized experient object
#' @param seFilteredSpliced seFilteredSpliced
#' @param exonsByReadClass exonsByReadClass
#' @param annotationGrangesList annotationGrangesList
#' @param filterSet1 filterSet1
#' @param min.exonOverlap min.exonOverlap
#' @param verbose verbose
#' @noRd
extdannotateUnsplicedReads <- function(se, seFilteredSpliced, exonsByReadClass,
                                       annotationGrangesList, filterSet1, min.exonOverlap, verbose) {
    start.ptm <- proc.time()
    if (any(rowData(se)$confidenceType == "unsplicedNew" & filterSet1)) {
        seFilteredUnspliced <-
            se[rowData(se)$confidenceType == "unsplicedNew" & filterSet1, ]
        exonsByReadClassUnspliced <- GenomicRanges::GRanges(
            seqnames = rowData(seFilteredUnspliced)$chr,
            ranges = IRanges(start = rowData(seFilteredUnspliced)$start,
                             end = rowData(seFilteredUnspliced)$end),
            strand = rowData(seFilteredUnspliced)$strand)
        partitioning <- PartitioningByEnd(seq_along(exonsByReadClassUnspliced),
                                          names = NULL)
        exonsByReadClassUnspliced$exon_rank <-
            rep(1, length(exonsByReadClassUnspliced))
        exonsByReadClassUnspliced$exon_endRank <-
            rep(1, length(exonsByReadClassUnspliced))
        exonsByReadClassUnspliced <-
            relist(exonsByReadClassUnspliced, partitioning)
        GenomeInfoDb::seqlevels(exonsByReadClassUnspliced) <- 
            unique(c(GenomeInfoDb::seqlevels(exonsByReadClassUnspliced),
                     GenomeInfoDb::seqlevels(annotationGrangesList)))
        mcols(seFilteredUnspliced)$GENEID <- NA
        mcols(seFilteredUnspliced)$readClassType <- "unsplicedNew"
        ## here: add filter to remove unspliced transcripts which overlap
        ## with known transcripts/high quality spliced transcripts
        overlapUnspliced <- findOverlaps(exonsByReadClassUnspliced,
                                         annotationGrangesList, minoverlap = min.exonOverlap,
                                         select = "first")
        seFilteredUnspliced <- seFilteredUnspliced[is.na(overlapUnspliced)]
        exonsByReadClassUnspliced <-
            exonsByReadClassUnspliced[is.na(overlapUnspliced)]
        ## combined spliced and unspliced Tx candidates
        seCombined <-
            SummarizedExperiment::rbind(seFilteredSpliced, seFilteredUnspliced)
        exonRangesCombined <- c(exonsByReadClass, exonsByReadClassUnspliced)
        names(exonRangesCombined) <- seq_along(exonRangesCombined)
    } else {
        seCombined <- seFilteredSpliced
        exonRangesCombined <- exonsByReadClass
        names(exonRangesCombined) <- seq_along(exonRangesCombined)
    }
    end.ptm <- proc.time()
    if (verbose) message("extended annotations for unspliced reads in ",
                         round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    return(list("seCombined" = seCombined, 
                "exonRangesCombined" = exonRangesCombined))
}

#' assigned read classes to annotated and new gene IDs
#' @param seCombined seCombined
#' @param exonRangesCombined exonRangesCombined
#' @param annotationGrangesList annotationGrangesList
#' @param min.exonOverlap min.exonOverlap
#' @param prefix prefix
#' @param verbose verbose
#' @noRd
assignGeneIDexonMatch <- function(seCombined, exonRangesCombined,
    annotationGrangesList,min.exonOverlap, prefix, verbose) {
    start.ptm <- proc.time()
    exonMatchGene <- findOverlaps(exonRangesCombined, annotationGrangesList,
        select = "arbitrary", minoverlap = min.exonOverlap)
    geneIdByExon <- rep(NA, length(exonRangesCombined))
    geneIdByExon[!is.na(exonMatchGene)] <- mcols(annotationGrangesList)$GENEID[
        exonMatchGene[!is.na(exonMatchGene)]]
    geneIdByExon[!is.na(mcols(seCombined)$GENEID)] <-
        mcols(seCombined)$GENEID[!is.na(mcols(seCombined)$GENEID)]
    exonMatchGene <- findOverlaps(exonRangesCombined[is.na(geneIdByExon)],
        exonRangesCombined[!is.na(geneIdByExon)],
        select = "arbitrary",minoverlap = min.exonOverlap)
    while (any(!is.na(exonMatchGene))) {
        geneIdByExon[is.na(geneIdByExon)][!is.na(exonMatchGene)] <-
            geneIdByExon[!is.na(geneIdByExon)][
                exonMatchGene[!is.na(exonMatchGene)]]
        exonMatchGene <- findOverlaps(exonRangesCombined[is.na(geneIdByExon)],
            exonRangesCombined[!is.na(geneIdByExon)],
            select = "arbitrary", minoverlap = min.exonOverlap)
    }
    mcols(seCombined)$GENEID[is.na(mcols(seCombined)$GENEID)] <-
        geneIdByExon[is.na(mcols(seCombined)$GENEID)]
    if (any(is.na(mcols(seCombined)$GENEID))) {
        newGeneIds <-
        assignNewGeneIds(exonRangesCombined[is.na(mcols(seCombined)$GENEID)],
        prefix = prefix, minoverlap = 5, ignore.strand = FALSE)
        mcols(seCombined)$GENEID[as.integer(newGeneIds$readClassId)] <- 
            newGeneIds$geneId
    }
    end.ptm <- proc.time()
    if (verbose) message("assigned read classes to annotated and
        new gene IDs in ", round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    return(seCombined)
}


#' Assign New Gene with Gene Ids
#' @param exByTx exByTx
#' @param prefix prefix, defaults to empty
#' @param minoverlap defaults to 5
#' @param ignore.strand defaults to FALSE
#' @noRd
assignNewGeneIds <- function(exByTx, prefix = "", minoverlap = 5,
                             ignore.strand = FALSE) {
    if (is.null(names(exByTx))) names(exByTx) <- seq_along(exByTx)
    
    exonSelfOverlaps <- findOverlaps(exByTx, exByTx, select = "all",
                                     minoverlap = minoverlap, ignore.strand = ignore.strand)
    hitObject <- as_tibble(exonSelfOverlaps) %>% arrange(queryHits, subjectHits)
    candidateList <- hitObject %>% group_by(queryHits) %>%
        filter(queryHits <= min(subjectHits), queryHits != subjectHits) %>%
        ungroup()
    filteredOverlapList <- hitObject %>% filter(queryHits < subjectHits)
    rm(list = c("exonSelfOverlaps", "hitObject"))
    gc(verbose = FALSE)
    length_tmp <- 1
    # loop to include overlapping read classes which are not in order
    while (nrow(candidateList) > length_tmp) {
        length_tmp <- nrow(candidateList)
        temp <- includeOverlapReadClass(candidateList, filteredOverlapList)
        candidateList <- rbind(temp, candidateList)
        while (nrow(temp)) { ## annotate transcripts by new gene id
            temp <- includeOverlapReadClass(candidateList, filteredOverlapList)
            candidateList <- rbind(temp, candidateList)} ## second loop
        tst <- candidateList %>% group_by(subjectHits) %>%
            mutate(subjectCount = n()) %>% group_by(queryHits) %>%
            filter(max(subjectCount) > 1) %>% ungroup()
        temp2 <- inner_join(tst, tst, by = c("subjectHits" = "subjectHits")) %>%
            filter(queryHits.x != queryHits.y) %>%
            mutate(
                queryHits = if_else(queryHits.x > queryHits.y,
                                    queryHits.y, queryHits.x),
                subjectHits = if_else(queryHits.x > queryHits.y,
                                      queryHits.x, queryHits.y)) %>%
            dplyr::select(queryHits, subjectHits) %>%
            distinct()
        candidateList <- distinct(rbind(temp2, candidateList))
    }
    candidateList <- candidateList %>% filter(!queryHits %in% subjectHits) %>%
        arrange(queryHits, subjectHits)
    idToAdd <-
        (which(!(seq_along(exByTx) %in% unique(candidateList$subjectHits))))
    candidateList <- rbind(candidateList, tibble(
        queryHits = idToAdd, subjectHits = idToAdd
    )) %>% arrange(queryHits, subjectHits) %>%
        mutate(geneId = paste("gene", prefix, ".", queryHits, sep = "")) %>%
        dplyr::select(subjectHits, geneId)
    candidateList$readClassId <- names(exByTx)[candidateList$subjectHits]
    candidateList <- dplyr::select(candidateList, readClassId, geneId)
    return(candidateList)
}


#' calculate distance between first and last exon matches
#' @param candidateList candidateList
#' @param filteredOverlapList filteredOverlapList
#' @noRd
includeOverlapReadClass <- function(candidateList, filteredOverlapList) {
    temp <- left_join(candidateList, filteredOverlapList,
                      by = c("subjectHits" = "queryHits")) %>%
        group_by(queryHits) %>%
        filter(!subjectHits.y %in% subjectHits, !is.na(subjectHits.y)) %>%
        ungroup() %>%
        dplyr::select(queryHits, subjectHits.y) %>%
        distinct() %>%
        dplyr::rename(subjectHits = subjectHits.y)
    return(temp)
}


#' calculate minimum equivalent classes for extended annotations
#' @param seCombined seCombined
#' @param annotationGrangesList annotationGrangesList
#' @param exonRangesCombined exonRangesCombined
#' @param prefix prefix
#' @param min.readFractionByGene min.readFractionByGene
#' @param min.sampleNumber min.sampleNumber
#' @param remove.subsetTx remove.subsetTx
#' @param verbose verbose
#' @noRd
filterTranscripts <- function(seCombined, annotationGrangesList,
                              exonRangesCombined, prefix, min.readFractionByGene,
                              min.sampleNumber, remove.subsetTx, verbose) {
    start.ptm <- proc.time() # (1) based on transcript usage
    countsTBL <- as_tibble(assays(seCombined)$counts, .name_repair = 
                               "unique") %>% mutate(geneId = mcols(seCombined)$GENEID) %>%
        group_by(geneId) %>% mutate_at(vars(-geneId), .funs = sum) %>%
        ungroup() %>% dplyr::select(-geneId)
    relCounts <- assays(seCombined)$counts / countsTBL
    filterTxUsage <- rowSums(relCounts >= min.readFractionByGene, na.rm = TRUE
    ) >= min.sampleNumber
    seCombinedFiltered <- seCombined[filterTxUsage]
    exonRangesCombinedFiltered <- exonRangesCombined[filterTxUsage]
    if (remove.subsetTx) { # (2) based on compatiblity with annotations
        exonRangesCombinedFiltered <- exonRangesCombinedFiltered[!grepl(
            "compatible",mcols(seCombinedFiltered)$readClassType)]
        seCombinedFiltered <- seCombinedFiltered[!grepl("compatible",
                                                        mcols(seCombinedFiltered)$readClassType)]
    }# (3) remove transcripts with identical junctions to annotations
    extendedAnnotationRanges <- removeTranscriptsWIdenJunct(
        seCombinedFiltered, exonRangesCombinedFiltered, 
        annotationGrangesList, prefix)
    end.ptm <- proc.time()
    if (verbose) message("transcript filtering in ",
                         round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    start.ptm <- proc.time()
    geneListWithNewTx <- which(mcols(extendedAnnotationRanges)$GENEID %in%
                                   mcols(extendedAnnotationRanges)$GENEID[
                                       which(mcols(extendedAnnotationRanges)$newTxClass != "annotation")])
    minEqClasses <-
        getMinimumEqClassByTx(extendedAnnotationRanges[geneListWithNewTx])
    end.ptm <- proc.time()
    if (verbose) message("calculated minimum equivalent classes for
        extended annotations in ", round((end.ptm - start.ptm)[3] / 60, 1),
                         " mins.")
    mcols(extendedAnnotationRanges)$eqClass[geneListWithNewTx] <-
        minEqClasses$eqClass[match(names(extendedAnnotationRanges[
            geneListWithNewTx]), minEqClasses$queryTxId)]
    mcols(extendedAnnotationRanges) <- mcols(extendedAnnotationRanges)[, 
                                                                       c("TXNAME", "GENEID", "eqClass", "newTxClass")]
    return(extendedAnnotationRanges)
}



#' remove transcripts with identical junctions to annotations
#' @noRd
removeTranscriptsWIdenJunct <- function(seCombinedFiltered, 
    exonRangesCombinedFiltered, annotationGrangesList, prefix){
    exonRangesCombinedFiltered <- exonRangesCombinedFiltered[mcols(
        seCombinedFiltered)$readClassType != "equal"]
    seCombinedFiltered <-
        seCombinedFiltered[mcols(seCombinedFiltered)$readClassType != "equal"]
    # simplified classification, can be further improved for readibility
    mcols(seCombinedFiltered)$newTxClass <- 
        mcols(seCombinedFiltered)$readClassType
    mcols(seCombinedFiltered)$newTxClass[mcols(seCombinedFiltered)$readClassType
        == "unsplicedNew" & grepl("gene", mcols(seCombinedFiltered)$GENEID)] <-
        "newGene-unspliced"
    mcols(seCombinedFiltered)$newTxClass[mcols(seCombinedFiltered)$readClassType
        == "allNew" & grepl("gene", mcols(seCombinedFiltered)$GENEID)] <-
        "newGene-spliced"
    extendedAnnotationRanges <- exonRangesCombinedFiltered
    mcols(extendedAnnotationRanges) <-
        mcols(seCombinedFiltered)[, c("GENEID", "newTxClass")]
    if (length(extendedAnnotationRanges)) 
        mcols(extendedAnnotationRanges)$TXNAME <- paste0(
            "tx",prefix, ".", seq_along(extendedAnnotationRanges))
    names(extendedAnnotationRanges) <- mcols(extendedAnnotationRanges)$TXNAME
    annotationRangesToMerge <- annotationGrangesList
    mcols(annotationRangesToMerge)$newTxClass <- 
        rep("annotation",length(annotationRangesToMerge))
    extendedAnnotationRanges <-
        c(extendedAnnotationRanges, annotationRangesToMerge)
    return(extendedAnnotationRanges)
}




#' Estimate distance between read class and annotations
#' @param seReadClass seReadClass
#' @inheritParams bambu
#' @noRd
isore.estimateDistanceToAnnotations <- function(seReadClass,
    annotationGrangesList, min.exonDistance = 35,
    min.primarySecondaryDist = 5, min.primarySecondaryDistStartEnd = 100000, 
    additionalFiltering = FALSE, verbose = FALSE) {
    start.ptm <- proc.time()
    readClassTable <-
        as_tibble(rowData(seReadClass), rownames = "readClassId") %>%
        dplyr::select(readClassId, confidenceType)

    distTable <- calculateDistToAnnotation(rowRanges(seReadClass),
        annotationGrangesList, maxDist = min.exonDistance,
        primarySecondaryDist = min.primarySecondaryDist,
        primarySecondaryDistStartEnd = min.primarySecondaryDistStartEnd,
        ignore.strand = FALSE)
    distTable$readCount <- assays(seReadClass)$counts[distTable$readClassId, ] 

    if (additionalFiltering) 
        distTable <- left_join(distTable, dplyr::select(readClassTable,
            readClassId, confidenceType), by = "readClassId") %>%
            mutate(relativeReadCount = readCount / txNumberFiltered)
    distTable <- dplyr::select(distTable, annotationTxId, readClassId,
                                readCount, compatible, equal,dist)
    distTable <- left_join(distTable, as_tibble(mcols(annotationGrangesList)[,
            c("TXNAME", "GENEID")]),by = c("annotationTxId" = "TXNAME"))

    end.ptm <- proc.time()
    if (verbose) message("calculated distance table in ",
            round((end.ptm - start.ptm)[3] / 60, 1), " mins.")

    readClassTable <- addGeneIdsToReadClassTable(readClassTable, distTable, 
                        seReadClass, verbose)
    metadata(seReadClass) <- list(distTable = distTable)
    rowData(seReadClass) <- readClassTable
    return(seReadClass)
}

#' generate readClassTable
#' @noRd
addGeneIdsToReadClassTable <- function(readClassTable, distTable, 
                                       seReadClass, verbose){
    readClassTable$equal <- readClassTable$readClassId %in%
        unlist((filter(distTable, equal) %>%
                    dplyr::select(readClassId) %>% distinct()))
    readClassTable$compatible <- readClassTable$readClassId %in%
        unlist((filter(distTable, compatible) %>% distinct()))
    start.ptm <- proc.time()
    # assign read classes to genes based on the highest read count per gene
    readClassToGeneIdTable <- dplyr::select(distTable, readClassId, GENEID,
                                            readCount) %>% group_by(GENEID) %>%
        mutate(geneCount = sum(readCount)) %>% distinct() %>%
        group_by(readClassId) %>% filter(geneCount == max(geneCount)) %>%
        filter(row_number() == 1) %>% 
        dplyr::select(readClassId, geneId = GENEID) %>% ungroup()
    newGeneCandidates <- (!readClassTable$readClassId %in%
                              readClassToGeneIdTable$readClassId)
    readClassToGeneIdTableNew <-
        assignNewGeneIds(rowRanges(seReadClass)[newGeneCandidates],
                         prefix = ".unassigned", minoverlap = 5, ignore.strand = FALSE)
    readClassGeneTable <- rbind(readClassToGeneIdTable, 
                                readClassToGeneIdTableNew)
    readClassTable <- left_join(readClassTable, readClassGeneTable,
                                by = "readClassId") %>%
        dplyr::select(confidenceType, geneId, compatible, equal)
    end.ptm <- proc.time()
    if (verbose) message("added gene Ids for each read class ",
                         round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    return(readClassTable)
}
