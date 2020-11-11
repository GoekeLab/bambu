#' create transcript model for splice junction reads
#' @param readGrgList reads GRangesList
#' @param annotationGrangesList annotation GRangesList
#' @param unlisted_junctions unlisted_junctions
#' @param uniqueJunctions uniqueJunctions
#' @param stranded stranded
#' @param verbose verbose
#' @noRd
createModelforJunctionReads <- function(readGrgList, annotationGrangesList,
    unlisted_junctions, uniqueJunctions, stranded, verbose) {
    GenomeInfoDb::seqlevels(unlisted_junctions) <-
        GenomeInfoDb::seqlevels(readGrgList)
    GenomeInfoDb::seqlevels(uniqueJunctions) <- 
        GenomeInfoDb::seqlevels(readGrgList)
    uniqueAnnotatedIntrons <- unique(unlistIntrons(annotationGrangesList,
        use.names = FALSE, use.ids = FALSE))
    junctionTables <- junctionStrandCorrection(uniqueJunctions,
        unlisted_junctions, uniqueAnnotatedIntrons,
        stranded = stranded, verbose = verbose)
    uniqueJunctions <- junctionTables[[1]][, c("score", "spliceMotif",
        "spliceStrand", "junctionStartName", "junctionEndName",
        "startScore", "endScore", "id")]
    unlisted_junctions <- junctionTables[[2]]
    uniqueJunctions$annotatedJunction <- (!is.na(GenomicRanges::match(
        uniqueJunctions,uniqueAnnotatedIntrons)))
    uniqueJunctions$annotatedStart <- uniqueJunctions$junctionStartName %in%
        uniqueJunctions$junctionStartName[uniqueJunctions$annotatedJunction]
    uniqueJunctions$annotatedEnd <- uniqueJunctions$junctionEndName %in%
        uniqueJunctions$junctionEndName[uniqueJunctions$annotatedJunction]
    uniqueJunctions <- correctJunctionFromPrediction(uniqueJunctions, verbose)
    start.ptm <- proc.time()
    readClassListSpliced <- constructSplicedReadClassTables(
        uniqueJunctions = uniqueJunctions,
        unlisted_junctions = unlisted_junctions,
        readGrgList = readGrgList,
        stranded = stranded)
    end.ptm <- proc.time()
    if (verbose)
    message("Finished create transcript models (read classes) for reads with
    spliced junctions in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    return(readClassListSpliced)
}

#' find unique junctions
#' @noRd
findUniqueJunctions <- function(uniqueJunctions, junctionModel, verbose){
    predictSpliceSites <- predictSpliceJunctions(
        annotatedJunctions = uniqueJunctions,
        junctionModel = junctionModel,
        verbose = verbose)
    uniqueJunctions <- predictSpliceSites[[1]][, c(
        "score", "spliceMotif",
        "spliceStrand", "junctionStartName", "junctionEndName",
        "startScore", "endScore", "annotatedJunction",
        "annotatedStart", "annotatedEnd")]
    return(list("uniqueJunctions" = uniqueJunctions, 
                "predictSpliceSites" = predictSpliceSites))
}

#' correct junction from prediction
#' @param uniqueJunctions uniqueJunctions
#' @param verbose verbose
#' @noRd
correctJunctionFromPrediction <- function(uniqueJunctions, verbose) {
    start.ptm <- proc.time()
    if (sum(uniqueJunctions$annotatedJunction) > 5000 &
        sum(!uniqueJunctions$annotatedJunction) > 4000) {
        uniqJunctionsNmodels <-
            findUniqueJunctions(uniqueJunctions, NULL, verbose)
        uniqueJunctions <- uniqJunctionsNmodels$uniqueJunctions
        junctionModel <- uniqJunctionsNmodels$predictSpliceSites[[2]]
    } else {
        junctionModel <- standardJunctionModels_temp
        uniqJunctionsNmodels <-
            findUniqueJunctions(uniqueJunctions, junctionModel, verbose)
        uniqueJunctions <- uniqJunctionsNmodels$uniqueJunctions
        message("Junction correction with not enough data,
            precalculated model is used")
    }
    end.ptm <- proc.time()
    if (verbose) 
        message("Model to predict true splice sites built in ",
            round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    start.ptm <- proc.time()
    uniqueJunctions <- findHighConfidenceJunctions( junctions = uniqueJunctions,
        junctionModel = junctionModel, verbose = verbose)
    uniqueJunctions$mergedHighConfJunctionIdAll_noNA <- 
        uniqueJunctions$mergedHighConfJunctionId
    uniqueJunctions$mergedHighConfJunctionIdAll_noNA[
        is.na(uniqueJunctions$mergedHighConfJunctionId)] <- 
        names(uniqueJunctions[is.na(uniqueJunctions$mergedHighConfJunctionId)])
    uniqueJunctions$strand.mergedHighConfJunction <- 
        as.character(strand(
            uniqueJunctions[uniqueJunctions$mergedHighConfJunctionIdAll_noNA]))
    end.ptm <- proc.time()
    if (verbose) 
        message("Finished correcting junction based on set of high confidence
            junctions in ", round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    return(uniqueJunctions)
}

#' generate exonByReadClass
#' @noRd
generateExonsByReadClass <- function(readGrgList, annotationGrangesList, 
    unlisted_junctions, uniqueJunctions, stranded, verbose){
    GenomeInfoDb::seqlevels(readGrgList) <-
        unique(c(GenomeInfoDb::seqlevels(readGrgList),
        GenomeInfoDb::seqlevels(annotationGrangesList)))
    GenomeInfoDb::seqlevels(annotationGrangesList) <- 
        GenomeInfoDb::seqlevels(readGrgList)
    readClassListSpliced <- createModelforJunctionReads(
        readGrgList, annotationGrangesList, unlisted_junctions,
        uniqueJunctions, stranded, verbose)
    # seqlevels are made equal (added for chromosomes missing in any of them)
    start.ptm <- proc.time()
    singleExonReads <- unlist(readGrgList[elementNROWS(readGrgList) == 1],
        use.names = FALSE)
    mcols(singleExonReads)$id <- mcols(readGrgList[
        elementNROWS(readGrgList) == 1])$id
    referenceExons <- unique(c(GenomicRanges::granges(unlist(
        readClassListSpliced[mcols(readClassListSpliced)$confidenceType ==
            "highConfidenceJunctionReads" &
            mcols(readClassListSpliced)$strand.rc != "*"], use.names = FALSE)), 
            GenomicRanges::granges(unlist(annotationGrangesList,
            use.names = FALSE))))
    readClassListUnsplicedWithAnnotation <- constructUnsplicedReadClasses(
        granges = singleExonReads, grangesReference = referenceExons,
        confidenceType = "unsplicedWithin", stranded = stranded)
    singleExonReads <- singleExonReads[!mcols(singleExonReads)$id %in%
        readClassListUnsplicedWithAnnotation$readIds]
    referenceExons <- reduce(singleExonReads, ignore.strand = !stranded)
    readClassListUnsplicedReduced <- constructUnsplicedReadClasses(
        granges = singleExonReads, grangesReference = referenceExons,
        confidenceType = "unsplicedNew", stranded = stranded)
    end.ptm <- proc.time()
    if (verbose) message("Finished create single exon transcript models
        (read classes) in ", round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    exonsByReadClass <- c(readClassListSpliced,
        readClassListUnsplicedWithAnnotation$exonsByReadClass,
        readClassListUnsplicedReduced$exonsByReadClass)
    return(exonsByReadClass)
}

#' Isoform reconstruction using genomic alignments
#' @param readGrgList readGrgList
#' @param runName runName
#' @param stranded stranded
#' @param quickMode quickMode
#' @inheritParams bambu
#' @noRd
isore.constructReadClasses <- function(readGrgList,
    runName = "sample1", annotationGrangesList,
    genomeSequence = NULL, stranded = FALSE, ncore = 1, verbose = FALSE) {
    unlisted_junctions <- unlistIntrons(readGrgList,
        use.ids = TRUE, use.names = FALSE)
    start.ptm <- proc.time()
    uniqueJunctions <- createJunctionTable(unlisted_junctions,
        ncore = ncore, genomeSequence = genomeSequence)
    # all seqlevels should be consistent, and drop those not in uniqueJunctions
    if (!all(GenomeInfoDb::seqlevels(unlisted_junctions) %in% 
        GenomeInfoDb::seqlevels(uniqueJunctions))) {
        unlisted_junctions <- GenomeInfoDb::keepSeqlevels(unlisted_junctions,
            value = GenomeInfoDb::seqlevels(unlisted_junctions)[
                GenomeInfoDb::seqlevels(unlisted_junctions) %in%
            GenomeInfoDb::seqlevels(uniqueJunctions)], pruning.mode = "coarse")
        readGrgList <- GenomeInfoDb::keepSeqlevels(readGrgList,
            value = GenomeInfoDb::seqlevels(readGrgList)[ 
            GenomeInfoDb::seqlevels(readGrgList) %in%
            GenomeInfoDb::seqlevels(uniqueJunctions)], pruning.mode = "coarse")
    } # the seqleels will be made comparable for all ranges,
    # warning is shown if annotation is missing some
    if (!all(GenomeInfoDb::seqlevels(readGrgList) %in% 
        GenomeInfoDb::seqlevels(annotationGrangesList))) 
        message("not all chromosomes present in reference annotations,
            annotations might be incomplete. Please compare objects
            on the same reference")
    end.ptm <- proc.time()
    if (verbose) message("Finished creating junction list with splice motif
        in ", round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    exonsByReadClass <- generateExonsByReadClass(readGrgList,
        annotationGrangesList, unlisted_junctions, uniqueJunctions,
        stranded, verbose)
    counts <- matrix(mcols(exonsByReadClass)$readCount,
        dimnames = list(names(exonsByReadClass), runName))
    colDataDf <- DataFrame(name = runName, row.names = runName)
    mcols(exonsByReadClass) <- mcols(exonsByReadClass)[, c("chr.rc", 
        "strand.rc", "intronStarts", "intronEnds", "confidenceType")]
    se <- SummarizedExperiment(assays = SimpleList(counts = counts),
        rowRanges = exonsByReadClass, colData = colDataDf)
    return(se)
}

#' for creating counts and start matrices
#' @noRd
createReadMatrix <- function(value, rowData, readClassSe){
    readMatrix <- matrix(value,dimnames = list(seq_len(nrow(rowData)),
        rownames(colData(readClassSe))), ncol = nrow(colData(readClassSe)),
        nrow = nrow(rowData))
    return(readMatrix)
}

#' create SE object for spliced Tx
#' @param rowData.spliced rowData.spliced
#' @param readClassSeRef reference readClass SummarizedExperiment
#' @param readClassSe readClass SummarizedExperiment
#' @param colDataCombined colDataCombined
#' @noRd
createSEforSplicedTx <- function(rowData.spliced, readClassSeRef,
                                readClassSe, colDataCombined) {
    counts.splicedRef <- createReadMatrix(0, rowData.spliced, readClassSeRef)
    start.splicedRef <- createReadMatrix(NA, rowData.spliced, readClassSeRef)
    end.splicedRef <- start.splicedRef
    counts.splicedNew <- createReadMatrix(0, rowData.spliced, readClassSe)
    start.splicedNew <- createReadMatrix(NA, rowData.spliced, readClassSe)
    end.splicedNew <- start.splicedNew
    counts.splicedRef[!is.na(rowData.spliced$id.ref), ] <-
        as.matrix(assays(readClassSeRef)$counts[
            rowData.spliced$id.ref[!is.na(rowData.spliced$id.ref)],])
    start.splicedRef[!is.na(rowData.spliced$id.ref), ] <-
        as.matrix(assays(readClassSeRef)$start[
            rowData.spliced$id.ref[!is.na(rowData.spliced$id.ref)],])
    end.splicedRef[!is.na(rowData.spliced$id.ref), ] <-
        as.matrix(assays(readClassSeRef)$end[
            rowData.spliced$id.ref[!is.na(rowData.spliced$id.ref)],])
    counts.splicedNew[!is.na(rowData.spliced$id.new), ] <-
        as.matrix(assays(readClassSe)$counts[
            rowData.spliced$id.new[!is.na(rowData.spliced$id.new)],])
    start.splicedNew[!is.na(rowData.spliced$id.new), ] <-
        as.matrix(rowData.spliced[!is.na(rowData.spliced$id.new), "start.new"])
    end.splicedNew[!is.na(rowData.spliced$id.new), ] <-
        as.matrix(rowData.spliced[!is.na(rowData.spliced$id.new), "end.new"])
    counts.spliced <- cbind(counts.splicedRef, counts.splicedNew)
    start.spliced <- cbind(start.splicedRef, start.splicedNew)
    end.spliced <- cbind(end.splicedRef, end.splicedNew)
    rowData.spliced$start <- rowMins(start.spliced, na.rm = TRUE)
    rowData.spliced$end <- rowMaxs(end.spliced, na.rm = TRUE)
    rowData.spliced <- dplyr::select(rowData.spliced, chr, start,
        end, strand, intronStarts, intronEnds) %>%
        mutate(confidenceType = "highConfidenceJunctionReads")
    se.spliced <- SummarizedExperiment(
        assays = SimpleList(counts = counts.spliced, 
        start = start.spliced,end = end.spliced),
        rowData = rowData.spliced, colData = colDataCombined)
    return(se.spliced)
}

#' create SE object for spliced Tx
#' @param readClassSeRef reference readClass SummarizedExperiment
#' @param readClassSe readClass SummarizedExperiment
#' @param readClassSeTBL readClassSeTBL
#' @param unsplicedRangesRef unsplicedRangesRef
#' @param unsplicedRangesNew unsplicedRangesNew
#' @param combinedSingleExonRanges combinedSingleExonRanges
#' @param rowData.unspliced rowData.unspliced
#' @param stranded stranded
#' @param colDataCombined colDataCombined
#' @noRd
createSEforUnsplicedTx <- function(readClassSeRef, readClassSe,
    readClassSeTBL, unsplicedRangesRef,
    unsplicedRangesNew, combinedSingleExonRanges,
    colDataCombined, rowData.unspliced, stranded) {
    refTables <- prepSEforUnsplicedTx(unsplicedRangesRef,
        combinedSingleExonRanges, readClassSeRef, stranded)
    counts.unsplicedRefSum <- refTables$counts.unsplicedSum
    start.unsplicedRefSum <- refTables$start.unsplicedSum
    end.unsplicedRefSum <- refTables$end.unsplicedSum
    newTables <- prepSEforUnsplicedTx(unsplicedRangesNew, 
        combinedSingleExonRanges, readClassSe, stranded, readClassSeTBL)
    counts.unsplicedNewSum <- newTables$counts.unsplicedSum
    start.unsplicedNewSum <- newTables$start.unsplicedSum
    end.unsplicedNewSum <- newTables$end.unsplicedSum
    counts.unsplicedRef <-
        createReadMatrix(0, rowData.unspliced, readClassSeRef)
    start.unsplicedRef <-
        createReadMatrix(NA, rowData.unspliced, readClassSeRef)
    end.unsplicedRef <- start.unsplicedRef
    counts.unsplicedNew <- createReadMatrix(0, rowData.unspliced, readClassSe)
    start.unsplicedNew <- createReadMatrix(NA, rowData.unspliced, readClassSe)
    end.unsplicedNew <- start.unsplicedNew
    counts.unsplicedRef[counts.unsplicedRefSum$index, ] <-
        as.matrix(counts.unsplicedRefSum[, colnames(counts.unsplicedRef)])
    start.unsplicedRef[counts.unsplicedRefSum$index, ] <-
        as.matrix(start.unsplicedRefSum[, colnames(start.unsplicedRef)])
    end.unsplicedRef[counts.unsplicedRefSum$index, ] <-
        as.matrix(end.unsplicedRefSum[, colnames(end.unsplicedRef)])
    counts.unsplicedNew[counts.unsplicedNewSum$index, ] <-
        as.matrix(counts.unsplicedNewSum[, colnames(counts.unsplicedNew)])
    start.unsplicedNew[counts.unsplicedNewSum$index, ] <-
        as.matrix(start.unsplicedNewSum[, "start"])
    end.unsplicedNew[counts.unsplicedNewSum$index, ] <-
        as.matrix(end.unsplicedNewSum[, "end"])
    counts.unspliced <- cbind(counts.unsplicedRef, counts.unsplicedNew)
    start.unspliced <- cbind(start.unsplicedRef, start.unsplicedNew)
    start.unspliced[which(is.infinite(start.unspliced))] <- NA
    end.unspliced <- cbind(end.unsplicedRef, end.unsplicedNew)
    end.unspliced[which(is.infinite(end.unspliced))] <- NA
    se.unspliced <- SummarizedExperiment(
        assays = SimpleList(counts = counts.unspliced,
        start = start.unspliced, end = end.unspliced),
        rowData = rowData.unspliced, colData = colDataCombined)
    return(se.unspliced)
}

#' prepare SE for unspliced Tx
#' @param unsplicedRanges unsplicedRanges
#' @param combinedSingleExonRanges combinedSingleExonRanges
#' @param readClass readClass
#' @param stranded stranded
#' @param readClassSeTBL default NULL
#' @noRd
prepSEforUnsplicedTx <- function(unsplicedRanges,
    combinedSingleExonRanges, readClass,
    stranded, readClassSeTBL = NULL) {
    overlapToCombined <-
        findOverlaps(unsplicedRanges, combinedSingleExonRanges,
        type = "within", ignore.strand = !stranded, select = "first")
    counts.unsplicedSum <- as_tibble(assays(readClass)[["counts"]])[
        rowData(readClass)$confidenceType == "unsplicedNew",] %>%
        mutate(index = overlapToCombined) %>%
        group_by(index) %>%
        summarise_all(sum, na.rm = TRUE)
    if (is.null(readClassSeTBL)) {
        start.unsplicedSum <- as_tibble(assays(readClass)[["start"]])[
            rowData(readClass)$confidenceType == "unsplicedNew",] %>%
            mutate(index = overlapToCombined) %>%
            group_by(index) %>%
            summarise_all(min, na.rm = TRUE)
        end.unsplicedSum <- as_tibble(assays(readClass)[["end"]])[
            rowData(readClass)$confidenceType == "unsplicedNew",] %>%
            mutate(index = overlapToCombined) %>%
            group_by(index) %>%
            summarise_all(max, na.rm = TRUE)
    } else {
        start.unsplicedSum <- readClassSeTBL %>%
            filter(confidenceType == "unsplicedNew") %>%
            dplyr::select(start) %>%
            mutate(index = overlapToCombined) %>%
            group_by(index) %>%
            summarise_all(min, na.rm = TRUE)
        end.unsplicedSum <- readClassSeTBL %>%
            filter(confidenceType == "unsplicedNew") %>%
            dplyr::select(end) %>%
            mutate(index = overlapToCombined) %>%
            group_by(index) %>%
            summarise_all(max, na.rm = TRUE)
    }
    tableList <- list(
        "counts.unsplicedSum" = counts.unsplicedSum,
        "start.unsplicedSum" = start.unsplicedSum,
        "end.unsplicedSum" = end.unsplicedSum
    )
    return(tableList)
}

#' create ref from a readClassSe object if readClassSeRef is not provided
#' @noRd
createRefFromReadClassSE <- function(readClassSe){
    counts <- assays(readClassSe)$counts
    start <- matrix(min(start(rowRanges(readClassSe))),
                    dimnames = dimnames(counts))
    end <- matrix(max(end(rowRanges(readClassSe))),
                    dimnames = dimnames(counts))
    rowData <- as_tibble(rowData(readClassSe))
    rowData$start <- rowMins(start)
    rowData$end <- rowMaxs(end)
    rowData <- rowData %>% dplyr::select(chr = chr.rc, start, end,
                strand = strand.rc, intronStarts, intronEnds, confidenceType)
    readClassSeRef <- SummarizedExperiment(
        assays = SimpleList(counts = counts, start = start, end = end),
        rowData = rowData, colData = colData(readClassSe))
    return(readClassSeRef)
}

#' Combine transcript candidates across samples
#' @param readClassSe readClassSe
#' @param readClassRef readClassRef
#' @param stranded stranded
#' @param verbose verbose
#' @noRd
isore.combineTranscriptCandidates <- function(readClassSe,
    readClassSeRef = NULL, stranded = FALSE, verbose = FALSE) {
    if (is.null(readClassSeRef)) { # create ref from a readClassSe object
        readClassSeRef <- createRefFromReadClassSE(readClassSe)
        return(readClassSeRef)
    } else {
        colDataCombined <- rbind(colData(readClassSeRef), colData(readClassSe))
        readClassSeRefTBL <- as_tibble(rowData(readClassSeRef), rownames = "id")
        readClassSeTBL <- as_tibble(rowData(readClassSe), rownames = "id") %>%
            mutate(start = min(start(rowRanges(readClassSe))),
                    end = max(end(rowRanges(readClassSe))))
        rowData.spliced <- full_join(filter(readClassSeRefTBL,
            confidenceType == "highConfidenceJunctionReads"), filter(
            readClassSeTBL, confidenceType == "highConfidenceJunctionReads"),
            by = c("chr" = "chr.rc", "strand" = "strand.rc",
                "intronStarts", "intronEnds"), suffix = c(".ref", ".new"))
        # create SE objects for spliced and unspliced Tx
        se.spliced <- createSEforSplicedTx(rowData.spliced, readClassSeRef,
                                            readClassSe, colDataCombined)
        readClassSeRefTBL.unspliced <- 
            filter(readClassSeRefTBL,confidenceType == "unsplicedNew")
        readClassSeTBL.unspliced <-
            filter(readClassSeTBL,confidenceType == "unsplicedNew")
        unsplicedRangesRef <- GenomicRanges::GRanges(
            seqnames = readClassSeRefTBL.unspliced$chr,
            ranges = IRanges(start = readClassSeRefTBL.unspliced$start,
                            end = readClassSeRefTBL.unspliced$end),
            strand = readClassSeRefTBL.unspliced$strand)
        unsplicedRangesNew <- GenomicRanges::GRanges(
            seqnames = readClassSeTBL.unspliced$chr.rc,
            ranges = IRanges(start = readClassSeTBL.unspliced$start,
                            end = readClassSeTBL.unspliced$end),
            strand = readClassSeTBL.unspliced$strand.rc)
        combinedSingleExonRanges <- reduce(
            c(unsplicedRangesRef,unsplicedRangesNew), ignore.strand = !stranded)
        rowData.unspliced <- as_tibble(combinedSingleExonRanges) %>%
            mutate_if(is.factor, as.character) %>% dplyr::select(chr = seqnames,
            start, end, strand = strand) %>% mutate(intronStarts = NA, 
            intronEnds = NA, confidenceType = "unsplicedNew")
        se.unspliced <-
            createSEforUnsplicedTx(readClassSeRef,readClassSe, readClassSeTBL,
            unsplicedRangesRef, unsplicedRangesNew,combinedSingleExonRanges,
            colDataCombined, rowData.unspliced, stranded)
        se.combined <- SummarizedExperiment::rbind(se.spliced, se.unspliced)
        rownames(se.combined) <- seq_len(nrow(se.combined))
        return(se.combined)
    }
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

#' update classificationTable
#' @noRd
updateWIntronMatches <- function(unlistedIntrons, unlistedIntronsAnnotations,
    partitioning, classificationTable, annotationGrangesList,
    seFilteredSpliced, exonsByReadClass, min.exonDistance){
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
            seFilteredSpliced, annotationGrangesList, min.exonDistance)
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

#' extended annotations for spliced reads
#' @param exonsByReadClass exonsByReadClass
#' @param intronsByReadClass intronsByReadClass
#' @param annotationGrangesList annotationGrangesList
#' @param seFilteredSpliced seFilteredSpliced
#' @param min.exonDistance min.exonDistance
#' @param verbose verbose
#' @noRd
extdannotateSplicedReads <- function(exonsByReadClass,intronsByReadClass, 
    annotationGrangesList, seFilteredSpliced, min.exonDistance, verbose) {
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
        seFilteredSpliced, exonsByReadClass, min.exonDistance)
    mcols(seFilteredSpliced)$readClassType <-
        apply(classificationTable, 1, paste, collapse = "")
    return(seFilteredSpliced)
}

#' assign gene id by maximum match
#' @noRd
assignGeneIDbyMaxMatch <- function(unlistedIntrons,
    unlistedIntronsAnnotations, overlapsNewIntronsAnnotatedIntrons,
    exonsByReadClass, seFilteredSpliced, annotationGrangesList,
    min.exonDistance) {
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
                byReadClass$seFilteredSpliced, min.exonDistance, verbose
            )
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
                                readCount, compatible, equal)
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


## helper functions to correct junctions
#' calculate stranded read counts
#' @param uniqueJunctions uniqueJunctions
#' @param junctionMatchList junctionMatchList
#' @param genomeSequence genomeSequence
#' @param unstranded_unlisted_junctions unstranded_unlisted_junctions
#' @param unlisted_junction_granges unlisted_junction_granges
#' @noRd
calculateStrandedReadCounts <- function(uniqueJunctions,
    genomeSequence,unstranded_unlisted_junctions,
    unlisted_junction_granges) {
    junctionMatchList <- methods::as(findMatches(uniqueJunctions,
        unstranded_unlisted_junctions),"List")
    uniqueJunctions_score <- elementNROWS(junctionMatchList)
    junctionStrandList <- extractList(strand(unlisted_junction_granges),
        junctionMatchList)
    junctionSeqStart <- BSgenome::getSeq(genomeSequence,
        IRanges::shift(flank(uniqueJunctions,width = 2), 2))#shift from IRanges
    junctionSeqEnd <- BSgenome::getSeq(genomeSequence,
        IRanges::shift(flank(uniqueJunctions,width = 2, start = FALSE), -2))
    junctionMotif <- paste(junctionSeqStart, junctionSeqEnd, sep = "-")
    junctionStartName <- paste(seqnames(uniqueJunctions),start(uniqueJunctions),
        sep = ":")
    junctionEndName <- paste(seqnames(uniqueJunctions), end(uniqueJunctions),
        sep = ":")
    startScore <- as.integer(tapply(uniqueJunctions_score,
        junctionStartName, sum)[junctionStartName])
    endScore <- as.integer(tapply(uniqueJunctions_score,
        junctionEndName, sum)[junctionEndName])
    mcols(uniqueJunctions) <- DataFrame(
        score = uniqueJunctions_score,
        plus_score = sum(junctionStrandList == "+"),
        minus_score = sum(junctionStrandList == "-"),
        spliceMotif = junctionMotif,
        spliceStrand = spliceStrand(junctionMotif),
        junctionStartName = junctionStartName,
        junctionEndName = junctionEndName,
        startScore = startScore,
        endScore = endScore,
        id = seq_along(uniqueJunctions))
    strand(uniqueJunctions) <- uniqueJunctions$spliceStrand
    return(uniqueJunctions)
}
