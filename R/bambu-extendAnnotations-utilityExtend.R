#' Extend annotations
#' @inheritParams bambu
#' @noRd
isore.extendAnnotations <- function(combinedTranscripts, annotationGrangesList,
                                    remove.subsetTx = TRUE,
                                    min.sampleNumber = 1, NDR = NULL, min.exonDistance = 35, min.exonOverlap = 10,
                                    min.primarySecondaryDist = 5, min.primarySecondaryDistStartEnd = 5, 
                                    min.readFractionByEqClass = 0, fusionMode = FALSE,
                                    prefix = "Bambu", baselineFDR = 0.1, defaultModels = NULL, verbose = FALSE, 
                                    predictStart = FALSE, predictEnd = FALSE){
  combinedTranscripts <- filterTranscripts(combinedTranscripts, min.sampleNumber)
  if (nrow(combinedTranscripts) > 0) {
    group_var <- c("intronStarts","intronEnds","chr","strand","start","end",
                   "confidenceType","readCount", "maxTxScore", "maxTxScore.noFit", 
                   "maxIntronChainScore", "maxIntronChainScore.noFit", 
                   "maxTssScore", "maxTssScore.noFit", "maxTesScore", "maxTesScore.noFit", 
                   "firstExonGroup", "lastExonGroup", "tesId",
                   "startRegionId", "endRegionId", "compatible", "equal")
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
    } else{ 
        rowDataFilteredSpliced <- NULL
    }
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
      remove.subsetTx, min.readFractionByEqClass, baselineFDR, NDR, defaultModels, verbose, predictStart, predictEnd)
    message(paste0("Novel transcripts detected: ", sum(mcols(extendedAnnotationRanges)$novelTranscript)))
    message(paste0("Novel genes detected: ", length(unique(mcols(extendedAnnotationRanges)$GENEID[mcols(extendedAnnotationRanges)$novelGene]))))
    message(paste0("Low confidence transcripts excluded: ", length(metadata(extendedAnnotationRanges)$lowConfidenceTranscripts)))
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
  #combinedTranscripts = combinedTranscripts[filterSet,]
  combinedTranscripts$maxTxScore[!filterSet] = -1
  combinedTranscripts$maxTxScore.noFit[!filterSet] = -1
  return(combinedTranscripts)
}

#' calculate minimum equivalent classes for extended annotations
#' @importFrom dplyr select as_tibble %>% mutate_at mutate group_by 
#'     ungroup .funs .name_repair vars 
#' @noRd
filterTranscriptsByAnnotation <- function(rowDataCombined, annotationGrangesList,
                                          exonRangesCombined, prefix,  remove.subsetTx, 
                                          min.readFractionByEqClass, baselineFDR = 0.1, 
                                          NDR = NULL, defaultModels = NULL, verbose, 
                                          predictStart = FALSE, predictEnd = FALSE) {
  start.ptm <- proc.time() # (1) based on transcript usage
  
  #calculate relative read count before any filtering
  rowDataCombined <- group_by(rowDataCombined, GENEID) %>% mutate(relReadCount = readCount/sum(readCount))
  
  rowDataCombined$TXNAME[is.na(rowDataCombined$TXNAME)] <- paste0(
      prefix, "Tx", seq_len(sum(is.na(rowDataCombined$TXNAME))))
  names(exonRangesCombined) <- rowDataCombined$TXNAME

  #filter out subset reads
  if (remove.subsetTx) { # (1) based on compatibility with annotations
    notCompatibleIds <- (!grepl("compatible", rowDataCombined$readClassType) |
        rowDataCombined$readClassType == "equal:compatible") #keep equal for FDR calculation
    subsetTranscripts <- combindRowDataWithRanges(
        rowDataCombined[!notCompatibleIds,], 
        exonRangesCombined[!notCompatibleIds])
    rowDataCombined$maxTxScore[grepl("compatible", rowDataCombined$readClassType) &
        rowDataCombined$readClassType != "equal:compatible"] <- -1
    rowDataCombined$maxTxScore.noFit[grepl("compatible", rowDataCombined$readClassType) &
        rowDataCombined$readClassType != "equal:compatible"] <- -1
  }
  #(2) remove transcripts below NDR threshold/identical junctions to annotations
  rowDataCombined <- calculateNDROnTranscripts(rowDataCombined, 
                        useTxScore = length(annotationGrangesList)==0)

  if(length(annotationGrangesList)>0){ #only recommend an NDR if its possible to calculate an NDR
      NDR <- recommendNDR(rowDataCombined, baselineFDR, NDR, defaultModels, verbose)
  } else if(is.null(NDR)) {
          NDR <- 0.5
  }
  filterSet <- (rowDataCombined$NDR.tx <= NDR | rowDataCombined$readClassType == "equal:compatible")
  lowConfidenceTranscripts <- combindRowDataWithRanges(
        rowDataCombined[!filterSet,], 
        exonRangesCombined[!filterSet])
  exonRangesCombined <- exonRangesCombined[filterSet]
  rowDataCombined <- rowDataCombined[filterSet,]
 
  

  mcols(exonRangesCombined)$txid <- seq_along(exonRangesCombined)
  minEq <- getMinimumEqClassByTx(exonRangesCombined)$eqClassById
  rowDataCombined$relSubsetCount <- rowDataCombined$readCount/unlist(lapply(minEq, function(x){return(sum(rowDataCombined$readCount[x]))}))
  #post extend annotation filters applied here (currently only subset filter)
  if(min.readFractionByEqClass>0 & sum(filterSet)>0) { # filter out subset transcripts based on relative expression
    filterSet <- rowDataCombined$relSubsetCount > min.readFractionByEqClass
    exonRangesCombined <- exonRangesCombined[filterSet]
    rowDataCombined <- rowDataCombined[filterSet,]
  }
  if(sum(filterSet)==0 & length(annotationGrangesList)==0) stop(
    "WARNING - No annotations were provided. Please increase NDR threshold to use novel transcripts")
  if(sum(filterSet)==0) message("WARNING - No novel transcripts meet the given thresholds. Try a higher NDR.")
  # (3) combine novel transcripts with annotations
  extendedAnnotationRanges <- combindRowDataWithRanges(rowDataCombined, exonRangesCombined)
  extendedAnnotationRanges <- combineWithAnnotations(
    rowDataCombined, extendedAnnotationRanges, 
    annotationGrangesList, prefix, predictStart, predictEnd)
  minEqClasses <-
    getMinimumEqClassByTx(extendedAnnotationRanges) # get eqClasses
  if(!identical(names(extendedAnnotationRanges),minEqClasses$queryTxId)) warning('eq classes might be incorrect')
  mcols(extendedAnnotationRanges)$eqClassById <- minEqClasses$eqClassById
  extendedAnnotationRanges <- calculateRelSubsetCount(extendedAnnotationRanges, minEqClasses$eqClassById, min.readFractionByEqClass)
  mcols(extendedAnnotationRanges) <- mcols(extendedAnnotationRanges)[, 
                 c("TXNAME", "GENEID", "NDR.tx", "NDR.ic", "NDR.tss", "NDR.tes", "novelGene", "novelTranscript", 
                 "txClassDescription","readCount","relReadCount", 
                 "relSubsetCount", "txid", "eqClassById", "maxTxScore", "maxTxScore.noFit", 
                 "maxTssScore", "maxTssScore.noFit", "maxTesScore", "maxTesScore.noFit", 
                 "maxIntronChainScore", "maxIntronChainScore.noFit", "startRegionId", "endRegionId", "equalRc.Tss", "equalRc.Tes", "anno.Tss", "anno.Tes")]
  metadata(extendedAnnotationRanges)$NDRthreshold = NDR
  if (remove.subsetTx) metadata(extendedAnnotationRanges)$subsetTranscripts = subsetTranscripts
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
    combinedTranscripts <- combinedTranscripts[combinedTranscripts$maxTxScore.noFit >=0, ] #ignore filtered out read classes
    equal <- combinedTranscripts$readClassType == "equal:compatible"
    equal[is.na(equal)] <- FALSE
    #add envirnment so poly() works
    attr(defaultModels$lmNDR[["terms"]], ".Environment") <- new.env(parent = parent.env(globalenv()))
    baseline <- predict(defaultModels$lmNDR, newdata=data.frame(NDR=baselineFDR))
    attr(defaultModels$lmNDR[["terms"]], ".Environment") <- c()
    score <- combinedTranscripts$maxTxScore.noFit
    score[is.na(score)] <- 0
    NDR.txScores <- calculateNDR(score, equal)
    NDR.tx.rec <- predict(lm(NDR.txScores~poly(score,3,raw=TRUE)), newdata=data.frame(score=baseline))
    NDR.tx.rec <- round(NDR.tx.rec,3)
    if(NDR.tx.rec > 1) NDR.tx.rec <- 0.999
    if (NDR.tx.rec < 0) NDR.tx.rec <- 0
    if(verbose) message("Recommended NDR for baseline FDR of ", baselineFDR, " = ", NDR.tx.rec)
    if(NDR.tx.rec > 0.5){
        message("A high NDR threshold is being recommended by Bambu indicating high levels of novel transcripts, ",
        "limiting the performance of the trained model")
        message("We recommend training a new model on similiar but well annotated dataset if available ",
        "(https://github.com/GoekeLab/bambu/tree/master#Training-a-model-on-another-speciesdataset-and-applying-it)",
        ", or alternatively running Bambu with opt.discovery=list(fitReadClassModel=FALSE)")
    }
    
    #if users are using an NDR let them know if the recommended NDR is different
    if(is.null(NDR)) {
        NDR <- NDR.tx.rec
        message("Using a novel discovery rate (NDR) of: ", NDR)
    } else if(abs(NDR.tx.rec-NDR)>=0.1){
            message(paste0("For your combination of sample and reference annotations we recommend an NDR of ", NDR.tx.rec,
            ". You are currently using an NDR threshold of ", NDR, 
            ". A higher NDR is suited for samples where the reference annotations are poor and more novel transcripts are expected,", 
            "whereas a lower NDR is suited for samples with already high quality annotations"))
    }
    return(NDR)
}

recommendNDR.onAnnotations <- function(annotations, prefix = "Bambu", baselineFDR = 0.1, defaultModels2 = defaultModels2){
    mcols <- mcols(annotations)[!is.na(mcols(annotations)$maxTxScore),]
    equal <- !grepl(prefix, mcols$TXNAME)
    #add envirnment so poly() works
    attr(defaultModels2$lmNDR[["terms"]], ".Environment") <- new.env(parent = parent.env(globalenv()))
    baseline <- predict(defaultModels2$lmNDR, newdata=data.frame(NDR.tx=baselineFDR))
    attr(defaultModels2$lmNDR[["terms"]], ".Environment") <- c()
    score <- mcols$maxTxScore.noFit
    NDR.txScores <- calculateNDR(score, equal)
    NDR.tx.rec <- predict(lm(NDR.txScores~poly(score,3,raw=TRUE)), newdata=data.frame(score=baseline))
    NDR.tx.rec <- round(NDR.tx.rec,3)
    return(NDR.tx.rec)
}


#' Calculate NDR based on transcripts 
#' @noRd
calculateNDROnTranscripts <- function(combinedTranscripts, useTxScore = FALSE){
      # calculate and filter by NDR
    equal <- combinedTranscripts$readClassType == "equal:compatible"
    equal[is.na(equal)] <- FALSE
    if(sum(equal, na.rm = TRUE)<50 | sum(!equal, na.rm = TRUE)<50 | useTxScore){
          combinedTranscripts$NDR.tx <- 1 - combinedTranscripts$maxTxScore
          combinedTranscripts$NDR.ic <- 1 - combinedTranscripts$maxIntronChainScore
          combinedTranscripts$NDR.tss <- 1 - combinedTranscripts$maxTssScore
          combinedTranscripts$NDR.tes <- 1 - combinedTranscripts$maxTesScore
          if(!useTxScore) message("WARNING - Less than 50 TRUE or FALSE read classes ",
            "for NDR precision stabilization.")
          message("NDR will be approximated as: (1 - Transcript Model Prediction Score)")
    } else {
        combinedTranscripts$NDR.tx <- calculateNDR(combinedTranscripts$maxTxScore, equal)
        combinedTranscripts$NDR.ic <- calculateNDR(combinedTranscripts$maxIntronChainScore, equal)
        combinedTranscripts$NDR.tss <- calculateNDR(combinedTranscripts$maxTssScore, equal)
        combinedTranscripts$NDR.tes <- calculateNDR(combinedTranscripts$maxTesScore, equal)
    }
    
    combinedTranscripts$NDR.tx[combinedTranscripts$maxTxScore==-1] <- 1
    combinedTranscripts$NDR.ic[combinedTranscripts$maxIntronChainScore==-1] <- 1
    combinedTranscripts$NDR.tss[combinedTranscripts$maxTssScore==-1] <- 1
    combinedTranscripts$NDR.tes[combinedTranscripts$maxTesScore==-1] <- 1
    return(combinedTranscripts)
}

#' calculates the minimum NDR for each score 
#' @noRd
calculateNDR <- function(score, labels){
    scoreOrder <- order(score, decreasing = TRUE) 
    labels <- labels[scoreOrder]
    NDR <- cumsum(!labels)/(seq_len(length(labels))) #calculate NDR
    NDR <- rev(cummin(rev(NDR))) #flatten NDR so its never higher than a lower ranked RC
    return(NDR[order(scoreOrder)]) #return to original order
}

#' generate exon/intron ByReadClass objects
#' @importFrom GenomicRanges makeGRangesListFromFeatureFragments
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @importFrom dplyr select distinct 
#' @noRd
makeExonsIntronsSpliced <- function(transcriptsTibble,annotationSeqLevels){
  if(all(is.na(transcriptsTibble$intronStarts))){
      intronsByReadClass <- GRangesList()
  } else {    
      intronsByReadClass <- makeGRangesListFromFeatureFragments(
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
  
  
  #add filtering based one what we have in 
  classificationTable$equal[which(rowDataFilteredSpliced$equal == FALSE)] <- ""
  classificationTable$compatible[which(rowDataFilteredSpliced$compatible == 0)] <- ""
  
  
  
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
  
  # remove TXNAME, GENEID, equal and compatible for 
  idx_startEnd <- which(rowDataFilteredSpliced$firstExonGroup == 0 |
                          rowDataFilteredSpliced$lastExonGroup == 0)
  #if (length(idx_startEnd) > 0) {
  #  classificationTable$compatible[idx_startEnd] <- ""
  #  classificationTable$equal[idx_startEnd] <- ""
  #  rowDataFilteredSpliced$GENEID[idx_startEnd] <- NA
  #  rowDataFilteredSpliced$TXNAME[idx_startEnd] <- NA
  #}

  
  classificationTable$compatible[which(rowDataFilteredSpliced$compatible == 0)] = ""
  classificationTable$equal[which(rowDataFilteredSpliced$compatible == 0)] = ""
  rowDataFilteredSpliced$GENEID[which(rowDataFilteredSpliced$compatible == 0)] = NA
  rowDataFilteredSpliced$TXNAME[which(rowDataFilteredSpliced$compatible == 0)] = NA  

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
  #classificationTable <- updateWStartEnd(rowDataFilteredSpliced, classificationTable)

  classificationTable$compatible[which(rowDataFilteredSpliced$compatible == 0)] = ""
  classificationTable$equal[which(rowDataFilteredSpliced$compatible == 0)] = ""

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
      !distNewTxByQuery$firstExonMatch]] <- "newFirstExon"
    classificationTable$newFirstExon[
      classificationTable$newFirstJunction != "newFirstJunction"] <- ""
    classificationTable$newLastExon[distNewTxByQuery$queryHits[
      !distNewTxByQuery$lastExonMatch]] <- "newLastExon"
    classificationTable$newLastExon[
      classificationTable$newLastJunction != "newLastJunction"] <- ""
  }
  return(classificationTable)
}


#' update classificationTable by start and end
#' @importFrom GenomicRanges match
#' @noRd
updateWStartEnd <- function(rowDataSplicedTibble, classificationTable) {
  idx <- which(rowDataSplicedTibble$firstExonGroup == 0 | 
                 rowDataSplicedTibble$lastExonGroup == 0)
  if (length(idx) > 0) {
    classificationTable$compatible[idx] <- ""
    classificationTable$equal[idx] <- ""
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
      firstExonMatch = any(firstExonMatch),
      lastExonMatch = any(lastExonMatch),
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
  print("check point 1!!! after findSpliceOverlapsByDist")
  txToAnTableFiltered <- genFilteredAnTable(spliceOverlaps,
      primarySecondaryDist, primarySecondaryDistStartEnd, DistCalculated = FALSE)
  # (2) calculate splice overlap for any not in the list (new exon >= 35bp)
  setTMP <- unique(txToAnTableFiltered$queryHits)

  spliceOverlaps_rest <- findSpliceOverlapsByDist(exByTx[-setTMP],
                                                    exByTxRef, maxDist = 0, type = "any", firstLastSeparate = TRUE,
                                                    dropRangesByMinLength = FALSE, cutStartEnd = TRUE,
                                                    ignore.strand = ignore.strand)
  print("check point 2!!! after findSpliceOverlapsByDist")

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
        print("check point 3!!! after findSpliceOverlapsByDist")
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
               uniqueFirstExonLengthQuery + uniqueLastExonLengthQuery) %>%
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
      filter((uniqueFirstExonLengthQuery <= primarySecondaryDistStartEnd &
                uniqueLastExonLengthQuery <= primarySecondaryDistStartEnd) ==
               max(uniqueFirstExonLengthQuery <=
                     primarySecondaryDistStartEnd & uniqueLastExonLengthQuery <=
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
    ndrCols <- colnames(rowDataCombinedFiltered)[grepl("NDR", colnames(rowDataCombinedFiltered))]
    if(length(ndrCols) > 0){
    mcols(extendedAnnotationRanges) <-
      rowDataCombinedFiltered[, c("TXNAME", "GENEID", "novelGene", "novelTranscript", "txClassDescription","readCount", ndrCols, "startRegionId", "endRegionId",
                                  "maxTxScore", "maxTxScore.noFit", "maxTssScore", "maxTssScore.noFit", "maxTesScore", "maxTesScore.noFit", 
                                  "maxIntronChainScore", "maxIntronChainScore.noFit", "relReadCount")]
    } else{
    mcols(extendedAnnotationRanges) <-
      rowDataCombinedFiltered[, c("TXNAME", "GENEID", "novelGene", "novelTranscript", "txClassDescription","readCount", "startRegionId", "endRegionId",
                                  "maxTxScore", "maxTxScore.noFit", "maxTssScore", "maxTssScore.noFit", "maxTesScore", "maxTesScore.noFit", 
                                  "maxIntronChainScore", "maxIntronChainScore.noFit", "relReadCount")]
    }
    return(extendedAnnotationRanges)
}

#' combine annotations with predicted transcripts
#' @noRd
combineWithAnnotations <- function(rowDataCombinedFiltered, 
                                   extendedAnnotationRanges,annotationGrangesList, prefix, 
                                   predictStart = FALSE, predictEnd = FALSE){
  equalRanges <- rowDataCombinedFiltered[!(rowDataCombinedFiltered$novelTranscript),]
  #use 5' and 3' ends from data to overwrite annotation starts and ends
  mcols(annotationGrangesList)$equalRc.Tss <- NA
  mcols(annotationGrangesList)$equalRc.Tes <- NA
  mcols(annotationGrangesList)$anno.Tss <- NA
  mcols(annotationGrangesList)$anno.Tes <- NA
  
  equalRanges$equalRc.Tss <- ifelse(equalRanges$strand == "+",
                                    equalRanges$start,
                                    equalRanges$end)
  equalRanges$equalRc.Tes <- ifelse(equalRanges$strand == "+",
                                    equalRanges$end,
                                    equalRanges$start)
  equalRanges$anno.Tss <- start(getTss(selectStartExonsFromGrangesList(annotationGrangesList[equalRanges$TXNAME], exonNumber = 1), width = 1))
  equalRanges$anno.Tes <- end(getTes(selectEndExonsFromGrangesList(annotationGrangesList[equalRanges$TXNAME], exonNumber = 1), width = 1))
  #unique the equalRanges by TXNAME to avoid multiple updates for TSS and TES

  equalRanges_unique <- equalRanges %>%
    group_by(TXNAME) %>%
    summarise(
      equalRc.Tss = equalRc.Tss[which.max(readCount)],
      equalRc.Tes = equalRc.Tes[which.max(readCount)],
      anno.Tss = unique(anno.Tss),
      anno.Tes = unique(anno.Tes),
      maxTxScore = weightedMean(maxTxScore, readCount),
      maxTxScore.noFit = weightedMean(maxTxScore.noFit, readCount),
      maxIntronChainScore = weightedMean(maxIntronChainScore, readCount),
      maxIntronChainScore.noFit = weightedMean(maxIntronChainScore.noFit, readCount),
      maxTssScore = weightedMean(maxTssScore, readCount),
      maxTssScore.noFit = weightedMean(maxTssScore.noFit, readCount),
      maxTesScore = weightedMean(maxTesScore, readCount),
      maxTesScore.noFit = weightedMean(maxTesScore.noFit, readCount),
      readCount = sum(readCount),
      relReadCount = sum(relReadCount),
      intronStarts = list(unique(intronStarts)),
      intronEnds = list(unique(intronEnds)),
      chr = list(unique(chr)),
      strand = list(unique(strand)),
      confidenceType = list(unique(confidenceType)),
      .groups = "drop")


  if(predictStart == TRUE){
    annotationGrangesList[equalRanges_unique$TXNAME] <- updateStartEnd(annotationGrangesList[equalRanges_unique$TXNAME], 
                                                                equalRanges_unique$equalRc.Tss, 
                                                                feature = "tss")
  }
  
  if(predictEnd == TRUE){
    annotationGrangesList[equalRanges_unique$TXNAME] <- updateStartEnd(annotationGrangesList[equalRanges_unique$TXNAME], 
                                                                       equalRanges_unique$equalRc.Tes,
                                                                feature = "tes")
  }
  #remove extended ranges that are already present in annotation
  extendedAnnotationRanges <- extendedAnnotationRanges[rowDataCombinedFiltered$novelTranscript]
  annotationRangesToMerge <- annotationGrangesList
  if(length(annotationGrangesList)){
    mcols(annotationRangesToMerge)$readCount <- NA
    mcols(annotationRangesToMerge)$txClassDescription <- "annotation"
    mcols(annotationRangesToMerge)$novelTranscript <- FALSE
    mcols(annotationRangesToMerge)$novelGene <- FALSE
    mcols(annotationRangesToMerge)$NDR.tx <- NA
    mcols(annotationRangesToMerge)$NDR.ic <- NA
    mcols(annotationRangesToMerge)$NDR.tss <- NA
    mcols(annotationRangesToMerge)$NDR.tes <- NA
    mcols(annotationRangesToMerge)$maxTxScore <- NA
    mcols(annotationRangesToMerge)$maxTxScore.noFit <- NA
    mcols(annotationRangesToMerge)$maxTssScore <- NA
    mcols(annotationRangesToMerge)$maxTssScore.noFit <- NA
    mcols(annotationRangesToMerge)$maxIntronChainScore <- NA
    mcols(annotationRangesToMerge)$maxIntronChainScore.noFit <- NA
    mcols(annotationRangesToMerge)$startRegionId <- NA
    mcols(annotationRangesToMerge)$endRegionId <- NA
    mcols(annotationRangesToMerge)$equalRc.Tss <- NA
    mcols(annotationRangesToMerge)$equalRc.Tes <- NA
    mcols(annotationRangesToMerge)$anno.Tss <- NA
    mcols(annotationRangesToMerge)$anno.Tes <- NA
    
    mcols(extendedAnnotationRanges) <- mcols(extendedAnnotationRanges)[,colnames(mcols(extendedAnnotationRanges))]
    #copy over stats to annotations from read classes
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$NDR.tx <- equalRanges_unique$NDR.tx
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$NDR.ic <- equalRanges_unique$NDR.ic
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$NDR.tss <- equalRanges_unique$NDR.tss
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$NDR.tes <- equalRanges_unique$NDR.tes
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$readCount <- equalRanges_unique$readCount
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$relReadCount <- equalRanges_unique$relReadCount
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$maxTxScore <- equalRanges_unique$maxTxScore
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$maxTxScore.noFit <- equalRanges_unique$maxTxScore.noFit
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$maxIntronChainScore <- equalRanges_unique$maxIntronChainScore
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$maxIntronChainScore.noFit <- equalRanges_unique$maxIntronChainScore.noFit
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$maxTssScore <- equalRanges_unique$maxTssScore
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$maxTssScore.noFit <- equalRanges_unique$maxTssScore.noFit
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$startRegionId <- equalRanges_unique$startRegionId
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$endRegionId <- equalRanges_unique$endRegionId
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$equalRc.Tss <- equalRanges_unique$equalRc.Tss
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$equalRc.Tes <- equalRanges_unique$equalRc.Tes
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$anno.Tss <- equalRanges_unique$anno.Tss
    mcols(annotationRangesToMerge[equalRanges_unique$TXNAME])$anno.Tes <- equalRanges_unique$anno.Tes
  }
  if (length(extendedAnnotationRanges)) {
    mcols(extendedAnnotationRanges)$TXNAME <- rowDataCombinedFiltered[rowDataCombinedFiltered$novelTranscript,]$TXNAME
    names(extendedAnnotationRanges) <- mcols(extendedAnnotationRanges)$TXNAME
    extendedAnnotationRanges <-
      c(extendedAnnotationRanges, annotationRangesToMerge) # this will throw error in line 648-649 when extendedAnnotationRanges is empty 
    mcols(extendedAnnotationRanges)$txid <- seq_along(extendedAnnotationRanges)
  } else{
    extendedAnnotationRanges <- annotationRangesToMerge
    mcols(extendedAnnotationRanges)$txid <- seq_along(extendedAnnotationRanges)
    mcols(extendedAnnotationRanges)$relReadCount <- NA
  }
  return(extendedAnnotationRanges)
}

#' update the start or end of annotated tx based one the read class
updateStartEnd <- function(grlist, newSites, feature) {
  if (length(grlist) != length(newSites))
    stop("`grlist` and `newSites` must have the same length")  
  gr_unlisted <- unlist(grlist, use.names = FALSE)
  gr_partition <- PartitioningByEnd(grlist)
  first_idx <- start(gr_partition)
  last_idx  <- end(gr_partition) 
  strand_vec <- as.character(strand(gr_unlisted)[first_idx])
  if (feature == "tss") {
    # For + strand → update first exon start
    start(gr_unlisted)[first_idx[strand_vec == "+"]] <- newSites[strand_vec == "+"]
    # For - strand → update last exon end
    end(gr_unlisted)[last_idx[strand_vec == "-"]] <- newSites[strand_vec == "-"]
  }
  if (feature == "tes") {
    # For + strand → update last exon end
    end(gr_unlisted)[last_idx[strand_vec == "+"]] <- newSites[strand_vec == "+"]
    # For - strand → update first exon start
    start(gr_unlisted)[first_idx[strand_vec == "-"]] <- newSites[strand_vec == "-"]
  }
  new_grlist <- relist(gr_unlisted, gr_partition)
  mcols(new_grlist) <- mcols(grlist)
  return(new_grlist)
}

#' calculate relative subset read count after filtering (increase speed, subsets are not considered here)'
#' @noRd
calculateRelSubsetCount <- function(extendedAnnotationRanges, minEq, min.readFractionByEqClass){
  filter <- !is.na(mcols(extendedAnnotationRanges)$readCount)
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
#   if (additionalFiltering) 
#     distTable <- left_join(distTable, select(readClassTable,
#                                              readClassId, confidenceType), by = "readClassId") %>%
#     mutate(relativeReadCount = readCount / txNumberFiltered)

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

#' Function to change NDR threshold on extendedAnnotations
#' @title Function to change NDR threshold on extendedAnnotations
#' @description This function train a model for use on other data
#' @param extendedAnnotations A GRangesList object produced from bambu(quant = FALSE) or rowRanges(se)
#' @param NDR The maximum NDR for novel transcripts to be in extendedAnnotations (0-1). If not provided a recommended NDR is calculated.
#' @param includeRef A boolean which if TRUE will also filter out reference annotations based on their NDR
#' @param prefix A string which determines which transcripts are considered novel by bambu and will be filtered (by default = 'Bambu')
#' @param baselineFDR a value between 0-1. Bambu uses this FDR on the trained model to recommend an equivilent NDR threshold to be used for the sample. By default, a baseline FDR of 0.1 is used. This does not impact the analysis if an NDR is set.
#' @param defaultModels a bambu trained model object that bambu will use when fitReadClassModel==FALSE or the data is not suitable for training, defaults to the pretrained model in the bambu package
#' Output - returns a similiar GRangesList object with entries swapped into or out of metadata(extendedAnnotations)$lowConfidenceTranscripts
#' @details 
#' @return extendedAnnotations with a new NDR threshold
#' @export
setNDR <- function(extendedAnnotations, NDR = NULL, includeRef = FALSE, prefix = 'Bambu', baselineFDR = 0.1, defaultModels2 = defaultModels){
    #Check to see if the annotations/gtf are dervived from Bambu
    if(is.null(mcols(extendedAnnotations)$NDR)){
        warning("Annotations were not extended by Bambu (or the wrong prefix was provided). NDR can not be set")
        return(extendedAnnotations)
    }

    if(is.null(metadata(extendedAnnotations)$lowConfidenceTranscripts)) 
        metadata(extendedAnnotations)$lowConfidenceTranscripts = GRangesList()

    #recommend an NDR (needed when users read in Bambu GTF)
    if(is.null(NDR)){
        tempAnno <- c(extendedAnnotations, metadata(extendedAnnotations)$lowConfidenceTranscripts)
        NDR <- recommendNDR.onAnnotations(tempAnno, prefix = prefix, baselineFDR = baselineFDR, defaultModels2 = defaultModels2)
        message("Recommending a novel discovery rate (NDR) of: ", NDR)
    }

    #If reference annotations should be filtered too (note that reference annotations with no read support arn't filtered)
    if(includeRef){
        toRemove <- (!is.na(mcols(extendedAnnotations)$NDR) & mcols(extendedAnnotations)$NDR > NDR)
        toAdd <- !is.na(mcols(metadata(extendedAnnotations)$lowConfidenceTranscripts)$NDR) & 
            mcols(metadata(extendedAnnotations)$lowConfidenceTranscripts)$NDR <= NDR  
    } else {
        toRemove <- (mcols(extendedAnnotations)$NDR > NDR & 
            grepl(prefix, mcols(extendedAnnotations)$TXNAME))
        toAdd <- (mcols(metadata(extendedAnnotations)$lowConfidenceTranscripts)$NDR <= NDR & 
            grepl(prefix, mcols(metadata(extendedAnnotations)$lowConfidenceTranscripts)$TXNAME))     
    }
  
  temp <- c(metadata(extendedAnnotations)$lowConfidenceTranscripts[!toAdd], extendedAnnotations[toRemove])
  extendedAnnotations <- c(extendedAnnotations[!toRemove], metadata(extendedAnnotations)$lowConfidenceTranscripts[toAdd])
  metadata(extendedAnnotations)$lowConfidenceTranscripts <- temp

  mcols(extendedAnnotations)$txid <- seq_along(extendedAnnotations)
  minEqClasses <- getMinimumEqClassByTx(extendedAnnotations)
  mcols(extendedAnnotations)$eqClassById <- minEqClasses$eqClassById
  
  metadata(extendedAnnotations)$NDRthreshold <- NDR

  return(extendedAnnotations)
}


#' Extend annotations by clusters
#' @noRd
isore.extendAnnotations.clusters <- function(readClassList, annotations, clusters, NDR, isoreParameters, stranded, bpParameters, fusionMode, verbose = FALSE){
    message("--- Start extending annotations for clusters ---")
    #if clustering is a csv, create a list with the barcodes for each cluster
    #csv must have two cols with heading barcode, cluster
    if(!is.list(clusters)){
        clusters <- read.csv(clusters)
        clusters <- clusters %>% group_by(cluster) %>% summarise(barcodes = list(barcode))
        clusters <- clusters$cluster
        clusters <- clusters$barcodes
        names(clusters) <- clusters
    }
    annotations.clusters <- list()
    rcfs.clusters <- list()
    clusters.rc <- splitReadClassFilesByRC(readClassList[[1]])
    txScores <- c()
    for(i in seq_along(clusters)){
        ###TODO need to account for the sample name here which is added to the barcode
        index <- match(clusters[[i]],gsub('demultiplexed','',metadata(readClassList[[1]])$samples)) 
        index <- index[!is.na(index)]
        if(length(index)<20) next
        rcf.counts <- clusters.rc[,index]
        rcf.filt <- readClassList[[1]][rowSums(rcf.counts)>0,]
        rowData(rcf.filt)$readCount <- rowSums(rcf.counts)[rowSums(rcf.counts)>0]
        countsTBL <- calculateGeneProportion(counts=mcols(rcf.filt)$readCount,
                                            geneIds=mcols(rcf.filt)$GENEID)
        rowData(rcf.filt)$geneReadProp <- countsTBL$geneReadProp
        rowData(rcf.filt)$geneReadCount <- countsTBL$geneReadCount
        rowData(rcf.filt)$startSD <- 0
        rowData(rcf.filt)$endSD <- 0
        rowData(rcf.filt)$readCount.posStrand <- 0
        thresholdIndex <- which(rowData(rcf.filt)$readCount>=isoreParameters$min.readCount)
        model <- trainBambu(rcf.filt, verbose = verbose, min.readCount = isoreParameters$min.readCount)
        txScore <- getTranscriptScore(rowData(rcf.filt)[thresholdIndex,], model,
                                defaultModels)
        rowData(rcf.filt)$txScore <- rep(NA,nrow(rcf.filt))
        rowData(rcf.filt)$txScore[thresholdIndex] <- txScore
        #txScores = cbind(txScores, rowData(rcf.filt)$txScore)
        rcfs.clusters[[names(clusters)[i]]] <- rcf.filt
        annotations.clusters[[names(clusters)[i]]] <- bambu.extendAnnotations(list(rcf.filt), annotations, NDR,
                                isoreParameters, stranded, bpParameters, fusionMode, verbose)
    }
    if(length(rcfs.clusters)>0){
        print("--- Merging all individual clusters ---")
        annotations.clusters[["merged"]] <- bambu.extendAnnotations(rcfs.clusters, annotations, NDR,
            isoreParameters, stranded, bpParameters, fusionMode, verbose)
    }
      
    return(annotations.clusters)
}