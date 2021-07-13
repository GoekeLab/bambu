#' Extend annotations
#' @inheritParams bambu
#' @noRd
isore.extendAnnotations <- function(combinedTranscripts, annotationGrangesList,
                                    remove.subsetTx = TRUE,
                                    min.sampleNumber = 1, min.exonDistance = 35, min.exonOverlap = 10,
                                    min.primarySecondaryDist = 5, min.primarySecondaryDistStartEnd = 5, 
                                    prefix = "", verbose = FALSE){
  readIndex = combinedTranscripts$readIndex
  combinedTranscripts = combinedTranscripts$combinedTranscripts     
  combinedTranscripts$tempID = seq_len(nrow(combinedTranscripts))
  filterSet <- filterTranscriptsByRead(combinedTranscripts, min.sampleNumber)
  if (any(filterSet, na.rm = TRUE)) {
    transcriptsFiltered <- combinedTranscripts[filterSet,]
    group_var <- c("intronStarts","intronEnds","chr","strand","start","end",
                   "confidenceType","readCount", "tempID")
    rowDataTibble <- select(transcriptsFiltered,all_of(group_var))
    annotationSeqLevels <- seqlevels(annotationGrangesList)
    rowDataSplicedTibble <- filter(rowDataTibble,
                                   confidenceType == "highConfidenceJunctionReads")
    transcriptRanges <- makeExonsIntronsSpliced(
      rowDataSplicedTibble, annotationSeqLevels)
    confidenceTypeVec <- rowDataTibble$confidenceType
    rowDataFilteredSpliced <- addNewSplicedReadClasses(transcriptRanges,
                                                       rowDataSplicedTibble, annotationGrangesList, 
                                                       min.exonDistance, min.primarySecondaryDist,
                                                       min.primarySecondaryDistStartEnd, verbose)
    rowDataFilteredUnspliced <- rowDataTibble[which(confidenceTypeVec == "unsplicedNew"),]
    SEnRng <- addNewUnsplicedReadClasses(rowDataFilteredUnspliced, 
                                         rowDataFilteredSpliced, transcriptRanges$exons, 
                                         annotationGrangesList, min.exonOverlap, verbose)
    rowDataCombined <- SEnRng$rowDataCombined
    exonRangesCombined <- SEnRng$exonRangesCombined
    countsCombined <- SEnRng$countsCombined
    # assign gene IDs based on exon match
    rowDataCombined <- assignGeneIDexonMatch(
      rowDataCombined, exonRangesCombined, annotationGrangesList,
      min.exonOverlap, prefix, verbose)
    ## filter out transcripts
    extendedAnnotationRanges <- filterTranscriptsByAnnotation(
      rowDataCombined, annotationGrangesList, exonRangesCombined, prefix,
      remove.subsetTx, verbose)
    metadata(extendedAnnotationRanges) = list(readIndex=readIndex)
    return(extendedAnnotationRanges)
  } else {
    message("The current filtering criteria filters out all new read 
            classes, please consider less strigent critria!")
    return(annotationGrangesList)
  }
}

#' filter transcripts by read counts
#' @noRd
filterTranscriptsByRead <- function(combinedTranscripts, min.sampleNumber){
  if (nrow(combinedTranscripts) > 0) 
    filterSet <- combinedTranscripts$NSampleReadCount >= min.sampleNumber & (
      combinedTranscripts$NSampleReadProp >= min.sampleNumber) & (
        combinedTranscripts$NSampleGeneFDR >= min.sampleNumber) & (
          combinedTranscripts$NSampleTxFDR >= min.sampleNumber)
  # filter based on read count and transcript usage
  return(filterSet)
}

#' calculate minimum equivalent classes for extended annotations
#' @importFrom dplyr select as_tibble %>% mutate_at mutate group_by 
#'     ungroup .funs .name_repair vars 
#' @noRd
filterTranscriptsByAnnotation <- function(rowDataCombined, annotationGrangesList,
                                          exonRangesCombined, prefix,  remove.subsetTx, verbose) {
  start.ptm <- proc.time() # (1) based on transcript usage
  if (remove.subsetTx) { # (1) based on compatiblity with annotations
    notCompatibleIds <- which(!grepl("compatible",
                                     rowDataCombined$readClassType))
    exonRangesCombined <- exonRangesCombined[notCompatibleIds]
    rowDataCombined <- rowDataCombined[notCompatibleIds,]
  }# (2) remove transcripts with identical junctions to annotations
  extendedAnnotationRanges <- removeTranscriptsWIdenJunct(
    rowDataCombined, exonRangesCombined, 
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
                                                                     c("TXNAME", "GENEID", "eqClass", "newTxClass","readCount", "tempID")]
  return(extendedAnnotationRanges)
}



#' generate exon/intron ByReadClass objects
#' @importFrom GenomicRanges makeGRangesListFromFeatureFragments
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @importFrom dplyr select distinct 
#' @noRd
makeExonsIntronsSpliced <- function(transcriptsTibble,annotationSeqLevels){
  intronsByReadClass <- makeGRangesListFromFeatureFragments(
    seqnames = transcriptsTibble$chr,
    fragmentStarts = transcriptsTibble$intronStarts,
    fragmentEnds = transcriptsTibble$intronEnds,
    strand = transcriptsTibble$strand)
  names(intronsByReadClass) <- seq_along(intronsByReadClass)
  seqlevels(intronsByReadClass) <-
    unique(c(seqlevels(intronsByReadClass), 
             annotationSeqLevels))
  exonsByReadClass <- createExonByReadClass(
    transcriptsTibble, annotationSeqLevels)
  return(list("exons" = exonsByReadClass, 
              "introns" = intronsByReadClass))
}


#' create exonsByReadClass
#' @param seFilteredSpliced a SummarizedExperiment object
#' for filtered spliced reads
#' @param annotationGrangesList annotation GRangesList object
#' @importFrom GenomicRanges makeGRangesListFromFeatureFragments unlist narrow
#'     elementNROWS relist
#' @importFrom GenomeInfoDb seqlevels
#' @noRd
createExonByReadClass <- function(transcriptsTibble, annotationSeqLevels) {
  exonEndsShifted <- paste(transcriptsTibble$intronStarts,
                           as.integer(transcriptsTibble$end + 1), sep = ",")
  exonStartsShifted <- paste(as.integer(transcriptsTibble$start - 1),
                             transcriptsTibble$intronEnds, sep = ",")
  exonsByReadClass <- makeGRangesListFromFeatureFragments(
    seqnames = transcriptsTibble$chr,
    fragmentStarts = exonStartsShifted,
    fragmentEnds = exonEndsShifted,
    strand = transcriptsTibble$strand)
  exonsByReadClass <- narrow(exonsByReadClass, start = 2, end = -2)
  # correct junction to exon differences in coordinates
  names(exonsByReadClass) <- seq_along(exonsByReadClass)
  # add exon start and exon end rank
  unlistData <- unlist(exonsByReadClass, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(exonsByReadClass)),
                                    names = NULL)
  exon_rank <- lapply(width((partitioning)), seq, from = 1)
  # * assumes positive for exon ranking
  negative_strand <- which(transcriptsTibble$strand == "-")
  exon_rank[negative_strand] <- lapply(exon_rank[negative_strand], rev) 
  exon_endRank <- lapply(exon_rank, rev)
  unlistData$exon_rank <- unlist(exon_rank)
  unlistData$exon_endRank <- unlist(exon_endRank)
  exonsByReadClass <- relist(unlistData, partitioning)
  seqlevels(exonsByReadClass) <-
    unique(c(seqlevels(exonsByReadClass), annotationSeqLevels))
  return(exonsByReadClass)
}



#' extended annotations for spliced reads
#' @noRd
addNewSplicedReadClasses <- function(combinedTranscriptRanges, 
                                     rowDataFilteredSpliced, annotationGrangesList, min.exonDistance, 
                                     min.primarySecondaryDist, min.primarySecondaryDistStartEnd, verbose){
  start.ptm <- proc.time()
  exonsByReadClass <- combinedTranscriptRanges$exons
  intronsByReadClass <- combinedTranscriptRanges$introns
  ovExon <- 
    findSpliceOverlapsQuick(cutStartEndFromGrangesList(exonsByReadClass),
                            cutStartEndFromGrangesList(annotationGrangesList)) # slow
  classificationTable <- 
    data.frame(matrix("", nrow = length(exonsByReadClass), ncol = 9), 
               stringsAsFactors = FALSE)
  colnames(classificationTable) <- c("equal", "compatible", "newWithin",
                                     "newLastJunction","newFirstJunction","newJunction", "allNew", 
                                     "newFirstExon", "newLastExon")
  equalQhits <- queryHits(ovExon[mcols(ovExon)$equal])
  classificationTable$equal[equalQhits[!duplicated(equalQhits)]] <- "equal"
  compatibleQhits <- queryHits(ovExon[mcols(ovExon)$compatible])
  classificationTable$compatible[
    compatibleQhits[!duplicated(compatibleQhits)]] <- "compatible"
  classificationTable$compatible[classificationTable$equal == "equal"] <- ""
  # annotate with transcript and gene Ids
  compatibleShits <- subjectHits(ovExon[mcols(ovExon)$compatible])
  rowDataFilteredSpliced$GENEID <- NA
  rowDataFilteredSpliced$GENEID[
    compatibleQhits[!duplicated(compatibleQhits)]] <-
    mcols(annotationGrangesList)$GENEID[
      compatibleShits[!duplicated(compatibleQhits)]]
  # annotate with compatible gene id,
  equalShits <- subjectHits(ovExon[mcols(ovExon)$equal])
  rowDataFilteredSpliced$GENEID[equalQhits[!duplicated(equalQhits)]] <-
    mcols(annotationGrangesList[equalShits[!duplicated(equalQhits)]])$GENEID
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
                         rowDataFilteredSpliced, exonsByReadClass, min.exonDistance,
                         min.primarySecondaryDist, min.primarySecondaryDistStartEnd)
  rowDataFilteredSpliced$readClassType <-
    apply(classificationTable, 1, paste, collapse = "")
  end.ptm <- proc.time()
  if (verbose) message("extended annotations for spliced reads in ",
                       round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
  return(rowDataFilteredSpliced)
}


#' update classificationTable
#' @importFrom GenomicRanges match
#' @noRd
updateWIntronMatches <- function(unlistedIntrons, unlistedIntronsAnnotations,
                                 partitioning, classificationTable, annotationGrangesList,
                                 rowDataFilteredSpliced, exonsByReadClass, min.exonDistance,
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
  strandVec <- rowDataFilteredSpliced$strand
  intronMatchAny <- any(intronMatchesList)
  classificationTable$newLastJunction[which( strandVec == "+" & !lastJunctionMatch &
                                               intronMatchAny | (strandVec == "-" & (!firstJunctionMatch &
                                                                                       intronMatchAny)))] <- "newLastJunction"
  classificationTable$newFirstJunction[which(strandVec == "+" & !firstJunctionMatch &
                                               intronMatchAny | (strandVec == "-" & (!lastJunctionMatch &
                                                                                       intronMatchAny)))] <- "newFirstJunction"
  classificationTable$newJunction[(
    sum(!intronMatchesList) > !firstJunctionMatch + !lastJunctionMatch) &
      intronMatchAny] <- "newJunction"
  classificationTable$allNew[!intronMatchAny] <- "allNew"
  ## assign gene ids based on the max # of matching introns/splice junctions
  overlapsNewIntronsAnnotatedIntrons <- findOverlaps(unlistedIntrons,
                                                     unlistedIntronsAnnotations, type = "equal",
                                                     select = "all", ignore.strand = FALSE)
  if (length(overlapsNewIntronsAnnotatedIntrons)) {
    distNewTxByQuery <- assignGeneIDbyMaxMatch(
      unlistedIntrons,unlistedIntronsAnnotations,
      overlapsNewIntronsAnnotatedIntrons, exonsByReadClass,
      rowDataFilteredSpliced, annotationGrangesList, min.exonDistance,
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
#' @importFrom dplyr as_tibble %>% group_by summarise filter ungroup
#' @noRd
assignGeneIDbyMaxMatch <- function(unlistedIntrons,
                                   unlistedIntronsAnnotations, overlapsNewIntronsAnnotatedIntrons,
                                   exonsByReadClass, rowDataFilteredSpliced, annotationGrangesList,
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
  rowDataFilteredSpliced$GENEID[is.na(rowDataFilteredSpliced$GENEID)] <-
    geneIdByIntron[is.na(rowDataFilteredSpliced$GENEID)]
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
addNewUnsplicedReadClasses <- function(rowDataFilteredUnspliced,  rowDataFilteredSpliced,
                                       exonsByReadClass, annotationGrangesList, min.exonOverlap, verbose) {
  start.ptm <- proc.time()
  rowDataCombined <- rowDataFilteredSpliced
  exonRangesCombined <- exonsByReadClass
  names(exonRangesCombined) <- seq_along(exonRangesCombined)
  if (nrow(rowDataFilteredUnspliced)) {
    exonsByReadClassUnspliced <- GRanges(
      seqnames = rowDataFilteredUnspliced$chr,
      ranges = IRanges(start = rowDataFilteredUnspliced$start,
                       end = rowDataFilteredUnspliced$end),
      strand = rowDataFilteredUnspliced$strand)
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
    rowDataFilteredUnspliced$GENEID <- NA
    rowDataFilteredUnspliced$readClassType <- "unsplicedNew"
    ## add filter to remove unspliced transcripts which overlap
    ## with known transcripts/high quality spliced transcripts
    overlapUnspliced <- findOverlaps(exonsByReadClassUnspliced,
                                     annotationGrangesList, minoverlap = min.exonOverlap,
                                     select = "first")
    naOverlapUnspliced <- which(is.na(overlapUnspliced))
    if(length(naOverlapUnspliced)) {
      rowDataFilteredUnspliced <- 
        rowDataFilteredUnspliced[naOverlapUnspliced,]
      exonsByReadClassUnspliced <-
        exonsByReadClassUnspliced[naOverlapUnspliced]
      ## combined spliced and unspliced Tx candidates
      rowDataCombined <- 
        bind_rows(rowDataFilteredSpliced, rowDataFilteredUnspliced)
      exonRangesCombined <- c(exonsByReadClass, exonsByReadClassUnspliced)
      names(exonRangesCombined) <- seq_along(exonRangesCombined)
    }
  } 
  end.ptm <- proc.time()
  if (verbose) message("extended annotations for unspliced reads in ",
                       round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
  return(list("rowDataCombined" = rowDataCombined, 
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
assignGeneIDexonMatch <- function(rowDataCombined, exonRangesCombined,
                                  annotationGrangesList,min.exonOverlap, prefix, verbose) {
  start.ptm <- proc.time()
  exonMatchGene <- findOverlaps(exonRangesCombined, annotationGrangesList,
                                select = "arbitrary", minoverlap = min.exonOverlap)
  geneIdByExon <- rep(NA, length(exonRangesCombined))
  geneIdByExon[!is.na(exonMatchGene)] <- mcols(annotationGrangesList)$GENEID[
    exonMatchGene[!is.na(exonMatchGene)]]
  geneIdByExon[!is.na(rowDataCombined$GENEID)] <-
    rowDataCombined$GENEID[!is.na(rowDataCombined$GENEID)]
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
  isNaGeneIdBefore <- is.na(rowDataCombined$GENEID)
  rowDataCombined$GENEID[isNaGeneIdBefore] <-
    geneIdByExon[isNaGeneIdBefore]
  isNaGeneIdAfter <- is.na(rowDataCombined$GENEID)
  if (any(isNaGeneIdAfter)) {
    newGeneIds <- 
      assignNewGeneIds(exonRangesCombined[isNaGeneIdAfter],
                       prefix = prefix, minoverlap = 5, ignore.strand = FALSE)
    rowDataCombined$GENEID[as.integer(newGeneIds$readClassId)] <- 
      newGeneIds$geneId
  }
  end.ptm <- proc.time()
  if (verbose) message("assigned read classes to annotated and
        new gene IDs in ", round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
  return(rowDataCombined)
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
      mutate(queryHits = if_else(queryHits.x > queryHits.y,
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



#' remove transcripts with identical junctions to annotations
#' @noRd
removeTranscriptsWIdenJunct <- function(rowDataCombinedFiltered, 
                                        exonRangesCombinedFiltered,annotationGrangesList, prefix){
  # commented as notEqual is already filtered in the beginning
  # notEqual <- which(rowDataCombinedFiltered$readClassType != "equal")
  # exonRangesCombinedFiltered <- exonRangesCombinedFiltered[notEqual]
  # rowDataCombinedFiltered <- rowDataCombinedFiltered[notEqual,]
  # countsTibbleFiltered <- countsTibbleFiltered[notEqual,]
  # simplified classification, can be further improved for readibility
  rowDataCombinedFiltered$newTxClass <- rowDataCombinedFiltered$readClassType
  rowDataCombinedFiltered$newTxClass[rowDataCombinedFiltered$readClassType
                                     == "unsplicedNew" & grepl("gene", rowDataCombinedFiltered$GENEID)] <-
    "newGene-unspliced"
  rowDataCombinedFiltered$newTxClass[rowDataCombinedFiltered$readClassType
                                     == "allNew" & grepl("gene", rowDataCombinedFiltered$GENEID)] <-
    "newGene-spliced"
  extendedAnnotationRanges <- exonRangesCombinedFiltered
  mcols(extendedAnnotationRanges) <-
    rowDataCombinedFiltered[, c("GENEID", "newTxClass","readCount", "tempID")]
  if (length(extendedAnnotationRanges)) 
    mcols(extendedAnnotationRanges)$TXNAME <- paste0(
      "tx",prefix, ".", seq_along(extendedAnnotationRanges))
  names(extendedAnnotationRanges) <- mcols(extendedAnnotationRanges)$TXNAME
  annotationRangesToMerge <- annotationGrangesList
  mcols(annotationRangesToMerge)$readCount <- 
    rep(NA,length(annotationRangesToMerge))
  mcols(annotationRangesToMerge)$newTxClass <- 
    rep("annotation",length(annotationRangesToMerge))
  extendedAnnotationRanges <-
    c(extendedAnnotationRanges, annotationRangesToMerge)
  return(extendedAnnotationRanges)
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
