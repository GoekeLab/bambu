
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
    readClassSeRefTBL <- as_tibble(rowData(readClassSeRef))
    readClassSeTBL <- as_tibble(rowData(readClassSe)) %>%
      mutate(start = min(start(rowRanges(readClassSe))),
             end = max(end(rowRanges(readClassSe))))
    
     # readClassSeRefTBL <- as_tibble(rowData(readClassSeRef), rownames = "id")
     # readClassSeTBL <- as_tibble(rowData(readClassSe), rownames = "id") %>%
     #   mutate(start = min(start(rowRanges(readClassSe))),
     #          end = max(end(rowRanges(readClassSe))))
    rowData.spliced <- full_join(filter(readClassSeRefTBL,
                                        confidenceType == "highConfidenceJunctionReads"), 
                                 filter(readClassSeTBL, 
                                        confidenceType == "highConfidenceJunctionReads"),
                                 by = c("chr" = "chr.rc", "strand" = "strand.rc",
                                        "intronStarts", "intronEnds"), suffix = c(".ref", ".new"))
    # create SE objects for spliced and unspliced Tx
    se.spliced <- createSEforSplicedTx(rowData.spliced, readClassSeRef,
                                       readClassSe, colDataCombined)
    
    ## alternative using dplyr data tables
    #min.readCount = 2, 
    #min.readFractionByGene = 0.05, 
    #min.sampleNumber = 1
  #  readClassSeRefTBL <- as_tibble(rowData(readClassSeRef), counts= assays(readClassSeRef)$counts[,1])
    
    
    ##
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



#' create ref from a readClassSe object if readClassSeRef is not provided
#' @noRd
createRefFromReadClassSE <- function(readClassSe){
  rownames(readClassSe) <- NULL
  counts <- assays(readClassSe)$counts
  start <- matrix(min(start(rowRanges(readClassSe))),
                  dimnames = dimnames(counts))
  end <- matrix(max(end(rowRanges(readClassSe))),
                dimnames = dimnames(counts))
  rowData <- as_tibble(rowData(readClassSe)) %>%
    mutate(start = rowMins(start),
           end = rowMaxs(end)) %>%
    dplyr::select(chr = chr.rc, start, end, strand = strand.rc, intronStarts,
                  intronEnds, confidenceType, id=rcId)
  readClassSeRef <- SummarizedExperiment(assays = SimpleList(counts = counts,
                                                             start = start,
                                                             end = end),
                                         rowData = rowData, 
                                         colData = colData(readClassSe))
  return(readClassSeRef)
}

#' create ref from a readClassSe object if readClassSeRef is not provided
#' @noRd
createRefFromReadClassSEOri <- function(readClassSe){
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



#' for creating counts and start matrices
#' @noRd
createReadMatrix <- function(value, rowData, readClassSe){
  readMatrix <- matrix(value,dimnames = list(seq_len(nrow(rowData)),
                                             rownames(colData(readClassSe))), ncol = nrow(colData(readClassSe)),
                       nrow = nrow(rowData))
  return(readMatrix)
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
