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
    indices = readClassListSpliced$indices
    readIds = as.numeric(names(indices))
    readClassListSpliced = readClassListSpliced$readClassList
    readClasses = as.data.frame(mcols(readClassListSpliced))
    readMatrix=matrix(nrow=length(readGrgList),ncol=2)
    readMatrix[readIds,]=cbind(readClasses$readClassId[indices],indices)
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
    rcListUnsplAnno <- constructUnsplicedReadClasses(
        granges = singleExonReads, grangesReference = referenceExons,
        confidenceType = "unsplicedWithin", stranded = stranded)
    rcsUnsplAnno = rcListUnsplAnno$readClasses
    rcsUnsplAnno = rcsUnsplAnno[unique(names(rcsUnsplAnno))]
    indices = rcListUnsplAnno$indices
    readIds = as.numeric(names(indices))
    readMatrix[readIds,1] = names(unlist(rcsUnsplAnno))[indices]
    readMatrix[readIds,2]= indices + length(readClassListSpliced)
    singleExonReads <- singleExonReads[!mcols(singleExonReads)$id %in% readIds]
    referenceExons <- reduce(singleExonReads, ignore.strand = !stranded)
    rcListUnsplReduced <- constructUnsplicedReadClasses(
        granges = singleExonReads, grangesReference = referenceExons,
        confidenceType = "unsplicedNew", stranded = stranded)
    end.ptm <- proc.time()
    if (verbose) message("Finished create single exon transcript models
        (read classes) in ", round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    rcsUnsplReduced = rcListUnsplReduced$readClasses
    rcsUnsplReduced = rcsUnsplReduced[unique(names(rcsUnsplReduced))]
    exonsByReadClass <- c(readClassListSpliced, rcsUnsplAnno,
        rcsUnsplReduced)
    indices = rcListUnsplReduced$indices
    indices2 = indices + length(readClassListSpliced) + length(rcsUnsplAnno)
    readIds = as.numeric(names(indices))
    readMatrix[readIds,1]=names(unlist(rcsUnsplReduced))[indices]
    readMatrix[readIds,2]=indices2 
    mcols(readGrgList)$readClass = readMatrix[,1]
    mcols(readGrgList)$readClassIndex = readMatrix[,2]
    return(list(exonsByReadClass=exonsByReadClass, readGrgList = readGrgList))
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
    generateExonsByReadClassOutput <- generateExonsByReadClass(readGrgList,
        annotationGrangesList, unlisted_junctions, uniqueJunctions,
        stranded, verbose)
    exonsByReadClass <- generateExonsByReadClassOutput$exonsByReadClass
    readGrgList <- generateExonsByReadClassOutput$readGrgList
    rm(generateExonsByReadClassOutput)
    counts <- matrix(mcols(exonsByReadClass)$readCount,
        dimnames = list(names(exonsByReadClass), runName))
    colDataDf <- DataFrame(name = runName, row.names = runName)
    mcols(exonsByReadClass) <- mcols(exonsByReadClass)[, c("chr.rc", 
        "strand.rc", "start.rc", "end.rc", "intronStarts", 
        "intronEnds", "confidenceType")]
    se <- SummarizedExperiment(assays = SimpleList(counts = counts),
        rowRanges = exonsByReadClass, colData = colDataDf)
    return(list(se=se, readGrgList = readGrgList))
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
#' @importFrom dplyr select %>% mutate
#' @noRd
createSEforSplicedTx <- function(rowData.spliced, readClassSeRef,
                                readClassSe, colDataCombined) {
    counts.splicedRef <- createReadMatrix(0, rowData.spliced, readClassSeRef)
    start.splicedRef <- createReadMatrix(NA, rowData.spliced, readClassSeRef)
    end.splicedRef <- start.splicedRef
    counts.splicedNew <- createReadMatrix(0, rowData.spliced, readClassSe)
    start.splicedNew <- createReadMatrix(NA, rowData.spliced, readClassSe)
    end.splicedNew <- start.splicedNew
    idrefs <- which(!is.na(rowData.spliced$id.ref))
    rD_idrefs <- rowData.spliced$id.ref[idrefs]
    idnews <- which(!is.na(rowData.spliced$id.new))
    rD_idnews <- rowData.spliced$id.new[idnews]
    counts.splicedRef[idrefs, ] <- 
        as.matrix(assays(readClassSeRef)$counts[rD_idrefs,])
    start.splicedRef[idrefs, ] <-
        as.matrix(assays(readClassSeRef)$start[rD_idrefs,])
    end.splicedRef[idrefs, ] <-
        as.matrix(assays(readClassSeRef)$end[rD_idrefs,])
    counts.splicedNew[idnews, ] <-
        as.matrix(assays(readClassSe)$counts[rD_idnews,])
    start.splicedNew[idnews, ] <-
        as.matrix(rowData.spliced[idnews, "start.new"])
    end.splicedNew[idnews, ] <-
        as.matrix(rowData.spliced[idnews, "end.new"])
    counts.spliced <- cbind(counts.splicedRef, counts.splicedNew)
    start.spliced <- cbind(start.splicedRef, start.splicedNew)
    end.spliced <- cbind(end.splicedRef, end.splicedNew)
    strand_bias.spliced = getSplicedAssay("strand_bias")
    startSD.spliced = getSplicedAssay("startSD")
    endSD.spliced = getSplicedAssay("endSD")
    rowData.spliced$start <- rowMins(start.spliced, na.rm = TRUE)
    rowData.spliced$end <- rowMaxs(end.spliced, na.rm = TRUE)
    rowData.spliced <- dplyr::select(rowData.spliced, chr, start,
        end, strand, intronStarts, intronEnds) %>%
        mutate(confidenceType = "highConfidenceJunctionReads")
    se.spliced <- SummarizedExperiment(
        assays = SimpleList(counts = counts.spliced, 
        start = start.spliced,end = end.spliced,
        strand_bias=strand_bias.spliced,
        startSD=startSD.spliced, endSD=endSD.spliced),
        rowData = rowData.spliced, colData = colDataCombined)
    print("createSEforSplicedTx 2" )
    return(se.spliced)
}

#' create SE object for spliced Tx
#' @param assay character "assay name"
#' @noRd
getSplicedAssay = function(assay){
    counts.splicedRef <- matrix(0,
        dimnames = list(1:nrow(rowData.spliced),
        rownames(colData(readClassSeRef))),
        ncol = nrow(colData(readClassSeRef)),
        nrow = nrow(rowData.spliced))
    
    counts.splicedNew <- matrix(0,
        dimnames = list(1:nrow(rowData.spliced),
        rownames(colData(readClassSe))),
        ncol = nrow(colData(readClassSe)),
        nrow = nrow(rowData.spliced))

    counts.splicedRef[!is.na(rowData.spliced$id.ref), ] <- 
        as.matrix(assays(readClassSeRef)[[assay]][rowData.spliced$id.ref[
            !is.na(rowData.spliced$id.ref)],])
    counts.splicedNew[!is.na(rowData.spliced$id.new), ] <- 
        as.matrix(assays(readClassSe)[[assay]][rowData.spliced$id.new[
            !is.na(rowData.spliced$id.new)],])
    
    counts.spliced <- cbind(counts.splicedRef, counts.splicedNew)
    return(counts.spliced)
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
    refSum_index <- counts.unsplicedRefSum$index
    newSum_index <- counts.unsplicedNewSum$index
    counts.unsplicedRef[refSum_index, ] <-
        as.matrix(counts.unsplicedRefSum[, colnames(counts.unsplicedRef)])
    start.unsplicedRef[refSum_index, ] <-
        as.matrix(start.unsplicedRefSum[, colnames(start.unsplicedRef)])
    end.unsplicedRef[refSum_index, ] <-
        as.matrix(end.unsplicedRefSum[, colnames(end.unsplicedRef)])
    counts.unsplicedNew[newSum_index, ] <-
        as.matrix(counts.unsplicedNewSum[, colnames(counts.unsplicedNew)])
    start.unsplicedNew[newSum_index, ] <-
        as.matrix(start.unsplicedNewSum[, "start"])
    end.unsplicedNew[newSum_index, ] <-
        as.matrix(end.unsplicedNewSum[, "end"])
    counts.unspliced <- cbind(counts.unsplicedRef, counts.unsplicedNew)
    start.unspliced <- cbind(start.unsplicedRef, start.unsplicedNew)
    start.unspliced[which(is.infinite(start.unspliced))] <- NA
    end.unspliced <- cbind(end.unsplicedRef, end.unsplicedNew)
    end.unspliced[which(is.infinite(end.unspliced))] <- NA
    strand_bias.unspliced = getUnsplicedAssay("strand_bias", 
        counts.unsplicedRefSum$index, counts.unsplicedNewSum$index, sum)
    startSD.unspliced = getUnsplicedAssay("startSD", 
        counts.unsplicedRefSum$index, counts.unsplicedNewSum$index, mean)
    endSD.unspliced = getUnsplicedAssay("endSD", 
        counts.unsplicedRefSum$index, counts.unsplicedNewSum$index, mean)
    se.unspliced <- SummarizedExperiment(
        assays = SimpleList(counts = counts.unspliced,
        start = start.unspliced, end = end.unspliced,
        strand_bias = strand_bias.unspliced,
        startSD = startSD.unspliced, endSD = endSD.unspliced),
        rowData = rowData.unspliced, colData = colDataCombined)
    return(se.unspliced)
}

#helper function for other features during createSEforUnsplicedTx
getUnsplicedAssay = function(feature, index.unsplicedRefSum, 
    index.unsplicedNewSum, fun){
    feature.unsplicedRefSum <- 
    as_tibble(assays(readClassSeRef)[[feature]])[
        rowData(readClassSeRef)$confidenceType=='unsplicedNew',] %>%
    mutate(index=overlapRefToCombined) %>%
    group_by(index) %>%
    summarise_all(fun, na.rm=TRUE)
    unsplicedNew.index <- rowData(readClassSe)$confidenceType=='unsplicedNew'
    feature.unsplicedNewSum <- 
        as_tibble(assays(readClassSe)[[feature]])[unsplicedNew.index,] %>% 
    mutate(index=overlapNewToCombined) %>%
        group_by(index) %>%
        summarise_all(fun, na.rm=TRUE)
    feature.unsplicedRef <- matrix(0,
        dimnames = list(1:nrow(rowData.unspliced),
        rownames(colData(readClassSeRef))),
        ncol = nrow(colData(readClassSeRef)),
        nrow = nrow(rowData.unspliced))
    feature.unsplicedNew <- matrix(0,
        dimnames = list(1:nrow(rowData.unspliced),
        rownames(colData(readClassSe))),
        ncol = nrow(colData(readClassSe)),
        nrow = nrow(rowData.unspliced))
    feature.unsplicedRef[index.unsplicedRefSum, ] <- 
        as.matrix(feature.unsplicedRefSum[,colnames(feature.unsplicedRef)])
    feature.unsplicedNew[index.unsplicedNewSum, ] <- 
        as.matrix(feature.unsplicedNewSum[,colnames(feature.unsplicedNew)])

    feature.unspliced <- cbind(feature.unsplicedRef, feature.unsplicedNew)
    return(feature.unspliced)
}

#' prepare SE for unspliced Tx
#' @param unsplicedRanges unsplicedRanges
#' @param combinedSingleExonRanges combinedSingleExonRanges
#' @param readClass readClass
#' @param stranded stranded
#' @param readClassSeTBL default NULL
#' @importFrom dplyr as_tibble %>% mutate group_by summarise_all filter select
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
#' @importFrom dplyr %>% select
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
    strand_bias = assays(readClassSe)$strand_bias
    startSD = assays(readClassSe)$startSD
    endSD = assays(readClassSe)$endSD
    assaysList = SimpleList(counts=counts,
            start=start, end=end, strand_bias=strand_bias,
            startSD=startSD, endSD=endSD)
    readClassSeRef <- SummarizedExperiment(
        assays = assaysList,
        rowData = rowData, colData = colData(readClassSe))
    return(readClassSeRef)
}

#' Combine transcript candidates across samples
#' @param readClassSe readClassSe
#' @param readClassRef readClassRef
#' @param stranded stranded
#' @param verbose verbose
#' @importFrom GenomicRanges GRanges
#' @importFrom SummarizedExperiment rbind
#' @importFrom dplyr as_tibble mutate_if select mutate %>% filter
#' @noRd
isore.combineTranscriptCandidates <- function(readClassSe,
    readClassSeRef = NULL, stranded = FALSE, verbose = FALSE) {
    if (is.null(readClassSeRef)) { # create ref from a readClassSe object
        readClassSeRef <- createRefFromReadClassSE(readClassSe)
        return(readClassSeRef)
    } else {
        colDataCombined <- 
            rbind(colData(readClassSeRef), colData(readClassSe))
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
        unsplicedRangesRef <- GRanges(
            seqnames = readClassSeRefTBL.unspliced$chr,
            ranges = IRanges(start = readClassSeRefTBL.unspliced$start,
                            end = readClassSeRefTBL.unspliced$end),
            strand = readClassSeRefTBL.unspliced$strand)
        unsplicedRangesNew <- GRanges(
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
#' @importFrom GenomicRanges makeGRangesListFromFeatureFragments unlist narrow
#'     elementNROWS relist
#' @importFrom GenomeInfoDb seqlevels
#' @noRd
createExonByReadClass <- function(seFilteredSpliced, annotationGrangesList) {
    exonEndsShifted <- paste(rowData(seFilteredSpliced)$intronStarts,
        as.integer(rowData(seFilteredSpliced)$end + 1), sep = ",")
    exonStartsShifted <- paste(as.integer(rowData(seFilteredSpliced)$start - 1),
        rowData(seFilteredSpliced)$intronEnds, sep = ",")
    exonsByReadClass <- makeGRangesListFromFeatureFragments(
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
    seqlevels(exonsByReadClass) <-
        unique(c(seqlevels(exonsByReadClass), seqlevels(annotationGrangesList)))
    return(exonsByReadClass)
}

#' update classificationTable
#' @importFrom GenomicRanges match
#' @noRd
updateWIntronMatches <- function(unlistedIntrons, unlistedIntronsAnnotations,
    partitioning, classificationTable, annotationGrangesList,
    seFilteredSpliced, exonsByReadClass, min.exonDistance,
    min.primarySecondaryDist, min.primarySecondaryDistStartEnd){
    intronMatches <- match(
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

#' assign gene id by maximum match
#' @importFrom dplyr as_tibble %>% group_by summarise filter ungroup
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

#' extended annotations for unspliced reads
#' @param se a summarized experient object
#' @param seFilteredSpliced seFilteredSpliced
#' @param exonsByReadClass exonsByReadClass
#' @param annotationGrangesList annotationGrangesList
#' @param filterSet1 filterSet1
#' @param min.exonOverlap min.exonOverlap
#' @param verbose verbose
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @importFrom SummarizedExperiment rbind
#' @noRd
extdannotateUnsplicedReads <- function(se, seFilteredSpliced, exonsByReadClass,
    annotationGrangesList, filterSet1, min.exonOverlap, verbose) {
    start.ptm <- proc.time()
    if (any(rowData(se)$confidenceType == "unsplicedNew" & filterSet1)) {
        seFilteredUnspliced <-
            se[rowData(se)$confidenceType == "unsplicedNew" & filterSet1, ]
        exonsByReadClassUnspliced <- GRanges(
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
        seqlevels(exonsByReadClassUnspliced) <- 
            unique(c(seqlevels(exonsByReadClassUnspliced),
            seqlevels(annotationGrangesList)))
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
#' @importFrom dplyr select as_tibble %>% mutate_at mutate group_by 
#'     ungroup .funs .name_repair vars 
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
#' @importFrom GenomicRanges makeGRangesListFromFeatureFragments
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @noRd
exonsintronsByReadClass <- function(se, annotationGrangesList, filterSet1){
    seFilteredSpliced <- se[rowData(se)$confidenceType ==
        "highConfidenceJunctionReads" & filterSet1, ]
    mcols(seFilteredSpliced)$GENEID <- NA
    intronsByReadClass <- makeGRangesListFromFeatureFragments(
        seqnames = rowData(seFilteredSpliced)$chr,
        fragmentStarts = rowData(seFilteredSpliced)$intronStarts,
        fragmentEnds = rowData(seFilteredSpliced)$intronEnds,
        strand = rowData(seFilteredSpliced)$strand)
    names(intronsByReadClass) <- seq_along(intronsByReadClass)
    seqlevels(intronsByReadClass) <-
        unique(c(seqlevels(intronsByReadClass), 
        seqlevels(annotationGrangesList)))
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

#' generate readClassTable
#' @importFrom dplyr select filter distinct unlist group_by mutate %>% ungroup
#'     left_join
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
    readClassToGeneIdTable <- select(distTable, readClassId, GENEID,
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
#' @importFrom dplyr select left_join as_tibble mutate %>%
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
        distTable <- left_join(distTable, select(readClassTable,
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


#' Assign New Gene with Gene Ids
#' @param exByTx exByTx
#' @param prefix prefix, defaults to empty
#' @param minoverlap defaults to 5
#' @param ignore.strand defaults to FALSE
#' @importFrom dplyr as_tibble %>% filter select ungroup inner_join mutate
#'     distinct if_else arrange tibble
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
#' @importFrom dplyr select rename %>% left_join group_by filter 
#'     ungroup distinct
#' @noRd
includeOverlapReadClass <- function(candidateList, filteredOverlapList) {
    temp <- left_join(candidateList, filteredOverlapList,
                    by = c("subjectHits" = "queryHits")) %>%
        group_by(queryHits) %>%
        filter(!subjectHits.y %in% subjectHits, !is.na(subjectHits.y)) %>%
        ungroup() %>%
        dplyr::select(queryHits, subjectHits.y) %>%
        distinct() %>%
        rename(subjectHits = subjectHits.y)
    return(temp)
}


#' Calculate distance from read class to annotation
#' @param exByTx exByTx
#' @param exByTxRef exByTxRef
#' @param maxDist defaults to 35
#' @param primarySecondaryDist defaults to 5
#' @param ignore.strand defaults to FALSE
#' @importFrom dplyr ungroup %>%
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
#' @importFrom dplyr as_tibble group_by %>% mutate n arrange filter arrange
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