#' Extend annotations
#' @inheritParams bambu
#' @noRd
isore.extendAnnotations <- function(combinedTranscripts, annotationGrangesList,
                                    remove.subsetTx = TRUE,
                                    min.sampleNumber = 1, NDR = NULL, min.exonDistance = 35, min.exonOverlap = 10,
                                    min.primarySecondaryDist = 5, min.primarySecondaryDistStartEnd = 5, 
                                    min.readFractionByEqClass = 0, fusionMode = FALSE,
                                    prefix = "Bambu", baselineFDR = 0.1, defaultModels = NULL, verbose = FALSE){
  combinedTranscripts <- filterTranscripts(combinedTranscripts, min.sampleNumber)
  if (nrow(combinedTranscripts) > 0) {
    group_var <- c("intronStarts","intronEnds","chr","strand","start","end",
                   "confidenceType","readCount", "maxTxScore", "maxTxScore.noFit")
    rowDataTibble <- select(combinedTranscripts,all_of(group_var))
    annotationSeqLevels <- seqlevels(annotationGrangesList)
    rowDataSplicedTibble <- filter(rowDataTibble,
                                   confidenceType == "highConfidenceJunctionReads")
    transcriptRanges <- makeExonsIntronsSpliced(
      rowDataSplicedTibble, annotationSeqLevels)
    confidenceTypeVec <- rowDataTibble$confidenceType
    if(nrow(rowDataSplicedTibble)>0){
        rowDataFilteredSpliced <- addNewSplicedReadClasses(transcriptRanges,
                                                       rowDataSplicedTibble, annotationGrangesList, 
                                                       min.exonDistance, min.primarySecondaryDist,
                                                       min.primarySecondaryDistStartEnd, verbose)
    } else{ rowDataFilteredSpliced = NULL}
    rowDataFilteredUnspliced <- rowDataTibble[which(confidenceTypeVec == "unsplicedNew"),]
    SEnRng <- addNewUnsplicedReadClasses(rowDataFilteredUnspliced, 
                                         rowDataFilteredSpliced, transcriptRanges$exons, 
                                         annotationGrangesList, min.exonOverlap, verbose)
    rowDataCombined <- SEnRng$rowDataCombined
    exonRangesCombined <- SEnRng$exonRangesCombined
    # assign gene IDs based on exon match
    geneIds <- assignGeneIds(grl = exonRangesCombined,
                             annotations = annotationGrangesList,
                             min.exonOverlap = min.exonOverlap,
                             fusionMode = fusionMode)
    rowDataCombined$GENEID <- geneIds[,1]
    rowDataCombined$novelGene <- geneIds[,2]
    if(fusionMode) rowDataCombined$readClassType[geneIds[,3]] <- 'fusionTranscript'
    
    # ## filter out transcripts
    extendedAnnotationRanges <- filterTranscriptsByAnnotation(
      rowDataCombined, annotationGrangesList, exonRangesCombined, prefix,
      remove.subsetTx, min.readFractionByEqClass, baselineFDR, NDR, defaultModels, verbose)
    return(extendedAnnotationRanges)
  } else {
    message("The current filtering criteria filters out all new read 
            classes, please consider less strigent critria!")
    return(annotationGrangesList)
  }
}

#' filter transcripts by read counts
#' @noRd
filterTranscripts <- function(combinedTranscripts, min.sampleNumber){
  filterSet = NULL
  if (nrow(combinedTranscripts) > 0){
    #filter by read counts
    naTxscore <- is.na(combinedTranscripts$NSampleTxScore)
    if(any(naTxscore)) 
       combinedTranscripts$NSampleTxScore[naTxscore] <- 0
    filterSet <- combinedTranscripts$NSampleReadCount >= min.sampleNumber & (
        combinedTranscripts$NSampleTxScore >= min.sampleNumber) & (
        combinedTranscripts$NSampleReadProp >= min.sampleNumber)
  }
  combinedTranscripts = combinedTranscripts[filterSet,]
  return(combinedTranscripts)
}

#' calculate minimum equivalent classes for extended annotations
#' @importFrom dplyr select as_tibble %>% mutate_at mutate group_by 
#'     ungroup .funs .name_repair vars 
#' @noRd
filterTranscriptsByAnnotation <- function(rowDataCombined, annotationGrangesList,
                                          exonRangesCombined, prefix,  remove.subsetTx, 
                                          min.readFractionByEqClass, baselineFDR = 0.1, 
                                          NDR = NULL, defaultModels = NULL, verbose) {
  start.ptm <- proc.time() # (1) based on transcript usage
  
  #calculate relative read count before any filtering
  rowDataCombined <- group_by(rowDataCombined, GENEID) %>% mutate(relReadCount = readCount/sum(readCount))
  
  #filter out subset reads
  if (remove.subsetTx) { # (1) based on compatibility with annotations
    notCompatibleIds <- (!grepl("compatible", rowDataCombined$readClassType) |
        rowDataCombined$readClassType == "equal:compatible") #keep equal for FDR calculation
    # subsetTranscripts <- combindRowDataWithRanges(
    #     rowDataCombined[!notCompatibleIds,], 
    #     exonRangesCombined[!notCompatibleIds])
    exonRangesCombined <- exonRangesCombined[notCompatibleIds]
    rowDataCombined <- rowDataCombined[notCompatibleIds,]
  }
  #(2) remove transcripts below NDR threshold/identical junctions to annotations
  rowDataCombined = calculateNDROnTranscripts(rowDataCombined, 
                        useTxScore = length(annotationGrangesList)==0)
  if(length(annotationGrangesList)>0){ #only recommend an NDR if its possible to calculate an NDR
      NDR = recommendNDR(rowDataCombined, baselineFDR, NDR, defaultModels, verbose)
  } else {
      if(is.null(NDR)) NDR = 0.5
  }
  filterSet = (rowDataCombined$NDR <= NDR | rowDataCombined$readClassType == "equal:compatible")
  lowConfidenceTranscripts <- combindRowDataWithRanges(
        rowDataCombined[!filterSet,], 
        exonRangesCombined[!filterSet])
  exonRangesCombined <- exonRangesCombined[filterSet]
  rowDataCombined <- rowDataCombined[filterSet,]

  mcols(exonRangesCombined)$txid <- seq_along(exonRangesCombined)

  if(sum(filterSet==0) & length(annotationGrangesList)==0) stop(
    "WARNING - No annotations were provided. Please increase NDR threshold to use novel transcripts")
  if(sum(filterSet)==0) message("WARNING - No novel transcripts meet the given thresholds. Try a higher NDR.")
  # (3) combine novel transcripts with annotations
  extendedAnnotationRanges <- combindRowDataWithRanges(rowDataCombined, exonRangesCombined)
  extendedAnnotationRanges <- combineWithAnnotations(
    rowDataCombined, extendedAnnotationRanges, 
    annotationGrangesList, prefix)
  minEqClasses <-
    getMinimumEqClassByTx(extendedAnnotationRanges) # get eqClasses
  if(!identical(names(extendedAnnotationRanges),minEqClasses$queryTxId)) warning('eq classes might be incorrect')
  mcols(extendedAnnotationRanges)$eqClassById <- minEqClasses$eqClassById
  extendedAnnotationRanges <- calculateRelSubsetCount(extendedAnnotationRanges, minEqClasses$eqClassById, min.readFractionByEqClass)
  mcols(extendedAnnotationRanges) <- mcols(extendedAnnotationRanges)[, 
                 c("TXNAME", "GENEID", "NDR", "novelGene", "novelTranscript", "txClassDescription","readCount","relReadCount", "relSubsetCount", "txid", "eqClassById")]
  metadata(extendedAnnotationRanges)$NDRthreshold = NDR
  #if (remove.subsetTx) metadata(extendedAnnotationRanges)$subsetTranscripts = subsetTranscripts
  metadata(extendedAnnotationRanges)$lowConfidenceTranscripts = lowConfidenceTranscripts
  end.ptm <- proc.time()
  if (verbose) message("transcript filtering in ",
                       round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
  return(extendedAnnotationRanges)
}

#' calculates an expected NDR based on the annotations'
#' @noRd
recommendNDR <- function(combinedTranscripts, baselineFDR = 0.1, NDR = NULL, defaultModels = defaultModels, verbose = FALSE){
    if(verbose) message("-- Predicting annotation completeness to determine NDR threshold --")
    equal = combinedTranscripts$readClassType == "equal:compatible"
    equal[is.na(equal)] = FALSE

    #add envirnment so poly() works
    attr(defaultModels$lmNDR[["terms"]], ".Environment") <- new.env(parent = parent.env(globalenv()))
    baseline = predict(defaultModels$lmNDR, newdata=data.frame(NDR=baselineFDR))
    attr(defaultModels$lmNDR[["terms"]], ".Environment") = c()

    NDRscores = calculateNDR(combinedTranscripts$maxTxScore.noFit, equal)
    score = combinedTranscripts$maxTxScore.noFit
    NDR.rec = predict(lm(NDRscores~poly(score,3,raw=TRUE)), newdata=data.frame(score=baseline))
    NDR.rec = round(NDR.rec,3)
    if(verbose) message("Recommended NDR for baseline FDR of ", baselineFDR, " = ", NDR.rec)
    if (NDR.rec < 0) NDR.rec = 0
    if(NDR.rec > 0.5){
        message("A high NDR threshold is being recommended by Bambu indicating high levels of novel transcripts, ",
        "limiting the performance of the trained model")
        message("We recommend training a new model on similiar but well annotated dataset if available ",
        "(https://github.com/GoekeLab/bambu/tree/master#Training-a-model-on-another-speciesdataset-and-applying-it)",
        ", or alternatively running Bambu with opt.discovery=list(fitReadClassModel=FALSE)")
    }
    
    #if users are using an NDR let them know if the recommended NDR is different
    if(is.null(NDR)) 
    {
        NDR = NDR.rec
        message("Using a novel discovery rate (NDR) of: ", NDR)
    }
    else{
        if(abs(NDR.rec-NDR)>=0.1){
            message(paste0("For your combination of sample and reference annotations we recommend an NDR of ", NDR.rec,
            ". You are currently using an NDR threshold of ", NDR, 
            ". A higher NDR is suited for samples where the reference annotations are poor and more novel transcripts are expected,", 
            "whereas a lower NDR is suited for samples with already high quality annotations"))
        }
    }
    return(NDR)
}


#' Calculate NDR based on transcripts 
#' @noRd
calculateNDROnTranscripts <- function(combinedTranscripts, useTxScore = FALSE){
      # calculate and filter by NDR
    equal = combinedTranscripts$readClassType == "equal:compatible"
    equal[is.na(equal)] = FALSE
    if(sum(equal, na.rm = TRUE)<50 | sum(!equal, na.rm = TRUE)<50 | useTxScore){
          combinedTranscripts$NDR = 1 - combinedTranscripts$maxTxScore
          if(!useTxScore) message("WARNING - Less than 50 TRUE or FALSE read classes ",
            "for NDR precision stabilization.")
          message("NDR will be approximated as: (1 - Transcript Model Prediction Score)")
    } else combinedTranscripts$NDR = calculateNDR(combinedTranscripts$maxTxScore, equal)
    return(combinedTranscripts)
}

#' calculates the minimum NDR for each score 
#' @noRd
calculateNDR = function(score, labels){
    scoreOrder = order(score, decreasing = TRUE) 
    labels = labels[scoreOrder]
    NDR = cumsum(!labels)/(seq_len(length(labels))) #calculate NDR
    NDR = rev(cummin(rev(NDR))) #flatten NDR so its never higher than a lower ranked RC
    return(NDR[order(scoreOrder)]) #return to original order
}

#' generate exon/intron ByReadClass objects
#' @importFrom GenomicRanges makeGRangesListFromFeatureFragments
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @importFrom dplyr select distinct 
#' @noRd
makeExonsIntronsSpliced <- function(transcriptsTibble,annotationSeqLevels){
  if(all(is.na(transcriptsTibble$intronStarts))){
      intronsByReadClass = GRangesList()}
  else {    intronsByReadClass <- makeGRangesListFromFeatureFragments(
        seqnames = transcriptsTibble$chr,
        fragmentStarts = transcriptsTibble$intronStarts,
        fragmentEnds = transcriptsTibble$intronEnds,
        strand = transcriptsTibble$strand)
  }
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
  equalSubHits <- subjectHits(ovExon[mcols(ovExon)$equal])
  rowDataFilteredSpliced$TXNAME <- NA
  rowDataFilteredSpliced$TXNAME[
    equalQhits[!duplicated(equalQhits)]] <-
    mcols(annotationGrangesList)$TXNAME[
      equalSubHits[!duplicated(equalQhits)]]
  compatibleSubHits <- subjectHits(ovExon[mcols(ovExon)$compatible])
  rowDataFilteredSpliced$GENEID <- NA
  rowDataFilteredSpliced$GENEID[
    compatibleQhits[!duplicated(compatibleQhits)]] <-
    mcols(annotationGrangesList)$GENEID[
      compatibleSubHits[!duplicated(compatibleQhits)]]
  # annotate with compatible gene id,
  rowDataFilteredSpliced$GENEID[equalQhits[!duplicated(equalQhits)]] <-
    mcols(annotationGrangesList[equalSubHits[!duplicated(equalQhits)]])$GENEID
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
    apply(classificationTable, 1, function(x){paste(x[x!=""], collapse = ":")})
  rowDataFilteredSpliced$novelTranscript = TRUE
  rowDataFilteredSpliced$novelTranscript[classificationTable$equal=="equal"] = FALSE
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
      primarySecondaryDist, primarySecondaryDistStartEnd, DistCalculated = FALSE)
  # (2) calculate splice overlap for any not in the list (new exon >= 35bp)
  setTMP <- unique(txToAnTableFiltered$queryHits)

  spliceOverlaps_rest <- findSpliceOverlapsByDist(exByTx[-setTMP],
                                                    exByTxRef, maxDist = 0, type = "any", firstLastSeparate = TRUE,
                                                    dropRangesByMinLength = FALSE, cutStartEnd = TRUE,
                                                    ignore.strand = ignore.strand)
  if(length(spliceOverlaps_rest) > 0){
    txToAnTableRest <-
      genFilteredAnTable(spliceOverlaps_rest, primarySecondaryDist,primarySecondaryDistStartEnd,
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
                                  txToAnTableRest, txToAnTableRestStartEnd )
  } 
  txToAnTableFiltered <- txToAnTableFiltered %>% ungroup()
  txToAnTableFiltered$readClassId <-
    names(exByTx)[txToAnTableFiltered$queryHits]
  txToAnTableFiltered$annotationTxId <-
    names(exByTxRef)[txToAnTableFiltered$subjectHits]
  txToAnTableFiltered$txid <-
    mcols(exByTxRef)$txid[txToAnTableFiltered$subjectHits]
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
    rowDataFilteredUnspliced$TXNAME <- NA
    rowDataFilteredUnspliced$readClassType <- "unsplicedNew"
    rowDataFilteredUnspliced$novelGene <- TRUE
    rowDataFilteredUnspliced$novelTranscript <- TRUE
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

#' extract the important row range columns and add them to the ranges for final output
#' @noRd
combindRowDataWithRanges <- function(rowDataCombinedFiltered, exonRangesCombinedFiltered){
    rowDataCombinedFiltered$txClassDescription <- rowDataCombinedFiltered$readClassType
    rowDataCombinedFiltered$txClassDescription[rowDataCombinedFiltered$readClassType
                                       == "unsplicedNew" & rowDataCombinedFiltered$novelGene] <-
      "newGene-unspliced"
    rowDataCombinedFiltered$txClassDescription[rowDataCombinedFiltered$readClassType
                                       == "allNew" & rowDataCombinedFiltered$novelGene] <-
      "newGene-spliced"
    extendedAnnotationRanges <- exonRangesCombinedFiltered
    mcols(extendedAnnotationRanges) <-
      rowDataCombinedFiltered[, c("GENEID", "novelGene", "novelTranscript", "txClassDescription","readCount", "NDR",
                                  "relReadCount")]
    return(extendedAnnotationRanges)
}

#' combine annotations with predicted transcripts
#' @noRd
combineWithAnnotations <- function(rowDataCombinedFiltered, 
                                        extendedAnnotationRanges,annotationGrangesList, prefix){
    equalRanges = rowDataCombinedFiltered[!is.na(rowDataCombinedFiltered$TXNAME),]
    #remove extended ranges that are already present in annotation
    extendedAnnotationRanges <- extendedAnnotationRanges[is.na(rowDataCombinedFiltered$TXNAME)]
    annotationRangesToMerge <- annotationGrangesList
    if(length(annotationGrangesList)){
        mcols(annotationRangesToMerge)$readCount <- NA
        mcols(annotationRangesToMerge)$txClassDescription <- "annotation"
        mcols(annotationRangesToMerge)$novelTranscript <- FALSE
        mcols(annotationRangesToMerge)$novelGene <- FALSE
        mcols(annotationRangesToMerge)$NDR <- NA
        mcols(extendedAnnotationRanges) <- mcols(extendedAnnotationRanges)[,colnames(mcols(extendedAnnotationRanges))]
        #copy over stats to annotations from read classes
        mcols(annotationRangesToMerge[equalRanges$TXNAME])$NDR = equalRanges$NDR
        mcols(annotationRangesToMerge[equalRanges$TXNAME])$readCount = equalRanges$readCount
        mcols(annotationRangesToMerge[equalRanges$TXNAME])$relReadCount = equalRanges$relReadCount
        #mcols(annotationRangesToMerge[equalRanges$TXNAME])$relSubsetCount = equalRanges$relSubsetCount
    }
    if (length(extendedAnnotationRanges)) {
      mcols(extendedAnnotationRanges)$TXNAME <- paste0(
      prefix, "Tx", seq_along(extendedAnnotationRanges))
    names(extendedAnnotationRanges) <- mcols(extendedAnnotationRanges)$TXNAME
    extendedAnnotationRanges <-
      c(extendedAnnotationRanges, annotationRangesToMerge) # this will throw error in line 648-649 when extendedAnnotationRanges is empty 
    mcols(extendedAnnotationRanges)$txid <- seq_along(extendedAnnotationRanges)
    }else{
      extendedAnnotationRanges <- annotationRangesToMerge
      mcols(extendedAnnotationRanges)$txid <- seq_along(extendedAnnotationRanges)
      mcols(extendedAnnotationRanges)$relReadCount = NA
      #mcols(extendedAnnotationRanges)$relSubsetCount = NA
    }
  return(extendedAnnotationRanges)
}

#' calculate relative subset read count after filtering (increase speed, subsets are not considered here)'
#' @noRd
calculateRelSubsetCount = function(extendedAnnotationRanges, minEq, min.readFractionByEqClass){
  filter = !is.na(mcols(extendedAnnotationRanges)$readCount)
  mcols(extendedAnnotationRanges)$relSubsetCount <- NA
  mcols(extendedAnnotationRanges)$relSubsetCount[filter] <- 
    mcols(extendedAnnotationRanges)$readCount[filter]/
    unlist(lapply(minEq[filter], function(x){return(sum(mcols(extendedAnnotationRanges)$readCount[x], na.rm = TRUE))}))
  #post extend annotation filters applied here (currently only subset filter)
  if(min.readFractionByEqClass>0) { # filter out subset transcripts based on relative expression
    filterSet <- is.na(mcols(extendedAnnotationRanges)$relSubsetCount) | mcols(extendedAnnotationRanges)$relSubsetCount > min.readFractionByEqClass
    extendedAnnotationRanges <- extendedAnnotationRanges[filterSet]
  }
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
  distTable <- dplyr::select(distTable, annotationTxId, txid, readClassId,
      readCount, compatible, equal,dist)
  distTable <- left_join(distTable, as_tibble(mcols(annotationGrangesList)[, c("txid", "GENEID")]),
      by = c("txid" = "txid"))
  
  end.ptm <- proc.time()
  if (verbose) message("calculated distance table in ",
                       round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
  readClassTable <- addGeneIdsToReadClassTable(readClassTable, distTable, 
                                               seReadClass, verbose)
  distTable <- DataFrame(distTable)
  distTable$eqClassById <- as.list(mcols(annotationGrangesList)$eqClassById)[distTable$txid]
  
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
      assignGeneIdsNoReference(rowRanges(seReadClass)[newGeneCandidates], 
      prefix = "unassigned")
  readClassToGeneIdTableNew <- tibble(names(seReadClass)[newGeneCandidates], 
      readClassToGeneIdTableNew)
  colnames(readClassToGeneIdTableNew) <- c('readClassId', 'geneId')
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
